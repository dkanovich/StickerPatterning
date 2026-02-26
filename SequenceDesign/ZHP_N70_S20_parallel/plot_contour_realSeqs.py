#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Canonical + KDE in Z-space with TWO real-sequence overlays
(Real-sequence Z(H), Z(p) computed with EXACT logic from compare_two_files_Zp_vs_Zh_inset_fixed.py)
+ Canonical 5 sequences and user-defined EXTRA_SEQUENCES markers.

Outputs:
• ZH_vs_Zp_KDE_with_two_real_overlays_exact.svg
"""

import os, glob, json, pickle, pathlib, math, time
from pathlib import Path
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import seaborn as sns

# ─── CONFIG ──────────────────────────────────────────────────────────
RESULT_DIR    = "results"  # JSONs with precomputed Z(H), Z(p) used ONLY for KDE
H_LOOKUP_PATH = Path("../assign_sticks/H_Z_lookup_100k.pkl")
P_LOOKUP_PATH = Path("../assign_sticks/P_Z_lookup_100k.pkl")

# Real set A
FILE_REAL_A   = "../assign_sticks/binary_masks_nuclear.txt"
LABEL_A       = "Nuclear"
COLOR_A       = "C0"
MARKER_A      = "o"
OUTLINE_A     = True   # True → hollow with colored edge; False → filled

# Real set B (optional)
FILE_REAL_B   = "../assign_sticks/binary_masks_PLCD.txt"
LABEL_B       = "PLCD"
COLOR_B       = "C3"
MARKER_B      = "o"
OUTLINE_B     = False

# Canonical sequences (same 5 as your KDE script)
CANON_SEQS = {
    "uniform": "2212212221221222122122212212221221222122122212212221221222122122212212",
    "term":    "1111111111111111111122222222222222222222222222222222222222222222222222",
    "middle":  "2222222222222222222222222111111111111111111112222222222222222222222222",
    "extreme": "1111111111222222222222222222222222222222222222222222222222221111111111",
    "patch":   "2222211112222222222111122222222221111222222222211112222222222111122222",
}
CANON_COLORS  = dict(term="green", uniform="red", middle="blue", extreme="purple", patch="orange")
CANON_MARKERS = dict(uniform='o', term='o', middle='o', extreme='o', patch='o')

# User extras (same interface as your KDE script)
EXTRA_SEQUENCES = [
    {"name": "H2_P1",  "seq": "2211122222111222222222222222222112222221112222111222222111222211122222", "marker": "D", "color": "black"},
    {"name": "H0_P3",  "seq": "2211222211222221122222211222222112222221122222211222211222211222221122", "marker": "X", "color": "black"},
    {"name": "H-2_P5","seq": "1122221122221222211222211222211222212222122221122222122221122222112222", "marker": "*", "color": "black"},
    {"name": "4patch","seq": "2222222111112222222222221111122222222222211111222222222222111112222222", "marker": "^", "color": "black"},
]
# Optional file with lines: name<TAB or comma>sequence-of-'1'/'2'
EXTRA_SEQ_FILE = None

OUTFILE       = "ZH_vs_Zp_KDE_with_two_real_overlays_exact.svg"

# EXACT real-seq mapping
STICKER_CHARS = set("S1s")
SPACER_CHARS  = set("Xx02")

plt.rcParams.update({"font.family": "arial", "font.size": 16})

# ─── progress helpers (EXACT style) ──────────────────────────────────
def now(): return time.strftime("%H:%M:%S")
def progress(msg): print(f"[{now()}] {msg}", flush=True)

# ─── H metric (EXACT copy for REAL) ──────────────────────────────────
class StickerPattern:
    def __init__(self, seq: str): self.seq = seq

    @staticmethod
    def calc_H(counts):
        total = sum(counts)
        if total == 0: return 0.0
        p = np.array(counts, float) / total
        return float(np.sum(-1.0 * p * np.ma.log(p).filled(0)))

    def _H_one_shift(self, bL: int, off: int):
        L = len(self.seq); counts = []
        for j in range(0, L, bL):
            block = self.seq[j+off:j+off+bL]
            if len(block) < bL: block += self.seq[:bL-len(block)]
            counts.append(block.count("S"))
        return self.calc_H(counts) / math.log(len(counts))  # length-normalised

    def calc_H_ave(self):
        return float(np.mean([self._H_one_shift(4, o) for o in range(4)] +
                             [self._H_one_shift(5, o) for o in range(5)]))

# ─── patchiness p = E/(N-1) (EXACT) ──────────────────────────────────
def patchiness(mask: str) -> float:
    N = mask.count("S")
    if N <= 1: return 0.0
    E = sum(1 for i in range(len(mask)-1) if mask[i]=="S" and mask[i+1]=="S")
    return E / (N - 1)

# ─── utils (EXACT) ───────────────────────────────────────────────────
def to_mask(raw: str) -> str:
    out = []
    for ch in raw.strip():
        if ch in STICKER_CHARS: out.append("S")
        elif ch in SPACER_CHARS: out.append("X")
    return "".join(out)

def load_lookup(p: Path):
    if not p.exists(): raise FileNotFoundError(f"Missing lookup: {p}")
    return pickle.loads(p.read_bytes())

def z_from_lookup(val: float, L: int, k: int, lookup: dict, which: str):
    pair = lookup.get(L, {}).get(k)
    if pair is None: return None, f"no_lookup_{which}"
    mu, sig = pair
    if sig <= 0.0:   return None, f"sigma_zero_{which}"
    return (val - mu) / sig, None

def read_lines(path: str):
    with open(path, "r") as f: return [line.rstrip("\n") for line in f]

def compute_real_Z_from_txt(txt_path: str, lookup_H: dict, lookup_P: dict, label: str, progress_every=0):
    lines = read_lines(txt_path)
    kept = {"Z_H": [], "Z_P": []}
    ignored = {"empty":0, "no_symbols":0,
               "no_lookup_H":0, "sigma_zero_H":0,
               "no_lookup_P":0, "sigma_zero_P":0,
               "nan":0}
    total = len(lines)
    progress(f"{label}: {total} lines")
    for i, raw in enumerate(lines):
        if progress_every and i and (i % progress_every == 0):
            progress(f"{label}: {i}/{total} ({100*i/total:.1f}%)")
        s = raw.strip()
        if not s: ignored["empty"] += 1; continue
        mask = to_mask(s)
        L = len(mask)
        if L == 0: ignored["no_symbols"] += 1; continue
        k = mask.count("S")
        try:
            H = StickerPattern(mask).calc_H_ave()
            p = patchiness(mask)
        except Exception:
            ignored["nan"] += 1; continue
        zH, rH = z_from_lookup(H, L, k, lookup_H, "H")
        if zH is None: ignored[rH] += 1; continue
        zP, rP = z_from_lookup(p, L, k, lookup_P, "P")
        if zP is None: ignored[rP] += 1; continue
        if not (np.isfinite(zH) and np.isfinite(zP)):
            ignored["nan"] += 1; continue
        kept["Z_H"].append(float(zH)); kept["Z_P"].append(float(zP))
    for k in kept: kept[k] = np.array(kept[k], float)
    progress(f"{label}: kept={len(kept['Z_H'])}, ignored={sum(ignored.values())}/{total}  "
             f"(no_lookup_H={ignored['no_lookup_H']}, sigma_zero_H={ignored['sigma_zero_H']}, "
             f"no_lookup_P={ignored['no_lookup_P']}, sigma_zero_P={ignored['sigma_zero_P']}, "
             f"no_symbols={ignored['no_symbols']}, empty={ignored['empty']}, nan={ignored['nan']})")
    return kept["Z_H"], kept["Z_P"]

# ─── Canonical/Extras (compute Z from '1'/'2' strings like your KDE script) ───
def _calc_H_counts_12(counts):
    p_arr = np.array(counts, float) / sum(counts)
    return float(np.sum(-1 * p_arr * ma.log(p_arr).filled(0)))

def _H_once_12(seq: str, aa: str, bL: int, offset: int) -> float:
    N = len(seq); counts = []
    for j in range(0, N, bL):
        block = seq[j + offset : j + offset + bL]
        if len(block) < bL:
            block += seq[: bL - len(block)]
        counts.append(block.count(aa))
    return _calc_H_counts_12(counts) / np.log(len(counts))

def H_raw_12(seq: str, aa: str = '1') -> float:
    return float(np.mean([np.mean([_H_once_12(seq, aa, bL, off) for off in range(bL)]) for bL in (4, 5)]))

def patchiness_12(seq: str, aa: str = '1') -> float:
    N = seq.count(aa)
    if N <= 1: return 0.0
    adj = sum(1 for i in range(len(seq)-1) if seq[i]==aa and seq[i+1]==aa)
    return adj / (N - 1)

def _lookup_mu_sigma(lookup, L, k):
    if lookup is None: return None
    entry_L = lookup.get(L, lookup.get(str(L)))
    if not isinstance(entry_L, dict): return None
    entry_k = entry_L.get(k, entry_L.get(str(k)))
    if entry_k is None: return None
    mu, sig = float(entry_k[0]), max(float(entry_k[1]), 1e-12)
    return mu, sig

def _mc_baseline(L, k, aa, n_ref=5000, seed=123):
    rng = np.random.default_rng(seed)
    base = np.array([aa]*k + ['2']*(L-k), dtype='<U1')
    H_vals = np.empty(n_ref, dtype=np.float32)
    p_vals = np.empty(n_ref, dtype=np.float32)
    for i in range(n_ref):
        rng.shuffle(base)
        s = ''.join(base)
        H_vals[i] = H_raw_12(s, aa)
        p_vals[i] = patchiness_12(s, aa)
    mu_H, sig_H = float(H_vals.mean()), float(H_vals.std(ddof=0))
    mu_p, sig_p = float(p_vals.mean()), float(p_vals.std(ddof=0))
    return mu_H, max(sig_H,1e-12), mu_p, max(sig_p,1e-12)

def to_Z_from_raw(H_val, p_val, L, k, aa='1', LOOKUP_H=None, LOOKUP_P=None):
    mu_sig_H = _lookup_mu_sigma(LOOKUP_H, L, k)
    mu_sig_P = _lookup_mu_sigma(LOOKUP_P, L, k)
    if mu_sig_H is None or mu_sig_P is None:
        mu_H, sH, mu_p, sP = _mc_baseline(L, k, aa)
    else:
        mu_H, sH = mu_sig_H
        mu_p, sP = mu_sig_P
    return (H_val - mu_H)/sH, (p_val - mu_p)/sP

def compute_point_from_seq(seq: str, aa='1', LOOKUP_H=None, LOOKUP_P=None):
    L = len(seq); k = seq.count(aa)
    H = H_raw_12(seq, aa); p = patchiness_12(seq, aa)
    return to_Z_from_raw(H, p, L, k, aa, LOOKUP_H, LOOKUP_P)

def _try_load_extra_from_file(pathlike):
    out = []
    if not pathlike: return out
    p = pathlib.Path(pathlike)
    if not p.exists():
        print(f"[warn] EXTRA_SEQ_FILE not found: {p}")
        return out
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"): continue
        parts = [x.strip() for x in line.replace(",", "\t").split("\t") if x.strip()]
        if len(parts) >= 2:
            out.append({"name": parts[0], "seq": parts[1]})
    return out

# ─── Main ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Load lookups
    progress("Loading H and P lookups …")
    LOOKUP_H = load_lookup(H_LOOKUP_PATH)
    LOOKUP_P = load_lookup(P_LOOKUP_PATH)
    progress("Lookups loaded.")

    # Load solution JSONs for KDE (NO solution scatter)
    json_paths = glob.glob(os.path.join(RESULT_DIR, "*.json"))
    if not json_paths:
        raise FileNotFoundError(f"No JSON files found in {RESULT_DIR}/")
    Z_H_vals, Z_p_vals = [], []
    for path in json_paths:
        with open(path) as fh:
            data = json.load(fh)
        ZH = data.get("H_seq") or data.get("h_seq") or data.get("H")
        Zp = data.get("p_seq") or data.get("p")
        if ZH is None or Zp is None: continue
        try:
            Z_H_vals.append(float(ZH)); Z_p_vals.append(float(Zp))
        except Exception: continue
    Z_H_vals = np.asarray(Z_H_vals, float)
    Z_p_vals = np.asarray(Z_p_vals, float)
    progress(f"KDE source points: n={Z_H_vals.size}")

    # Compute real sets using EXACT original logic
    real_sets = []
    if FILE_REAL_A and Path(FILE_REAL_A).exists():
        ZH_A, Zp_A = compute_real_Z_from_txt(FILE_REAL_A, LOOKUP_H, LOOKUP_P, LABEL_A, progress_every=500)
        real_sets.append((ZH_A, Zp_A, LABEL_A, COLOR_A, MARKER_A, OUTLINE_A))
    else:
        progress(f"[warn] Real A not found or unset: {FILE_REAL_A}")
    if FILE_REAL_B and Path(FILE_REAL_B).exists():
        ZH_B, Zp_B = compute_real_Z_from_txt(FILE_REAL_B, LOOKUP_H, LOOKUP_P, LABEL_B, progress_every=500)
        real_sets.append((ZH_B, Zp_B, LABEL_B, COLOR_B, MARKER_B, OUTLINE_B))
    else:
        progress(f"[warn] Real B not found or unset: {FILE_REAL_B}")

    # Canonical and Extra points (display-only)
    canon_pts = []
    for name, seq in CANON_SEQS.items():
        ZH, Zp = compute_point_from_seq(seq, '1', LOOKUP_H, LOOKUP_P)
        canon_pts.append((name, ZH, Zp, CANON_MARKERS[name], CANON_COLORS[name]))

    EXTRA_SEQUENCES = EXTRA_SEQUENCES + _try_load_extra_from_file(EXTRA_SEQ_FILE)
    extra_pts = []
    for item in EXTRA_SEQUENCES:
        name = item.get("name", "extra")
        seq  = item.get("seq", "")
        if not seq or set(seq) - set("12"):
            print(f"[warn] Skipping invalid extra sequence '{name}': must be '1'/'2' only.")
            continue
        ZH, Zp = compute_point_from_seq(seq, '1', LOOKUP_H, LOOKUP_P)
        marker = item.get("marker", "*")
        color  = item.get("color", None)
        extra_pts.append((name, ZH, Zp, marker, color))

    # ─── Plot ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(6.9, 6.9))

    # KDE only (no solution scatter points)
    if Z_H_vals.size > 1:
        sns.kdeplot(x=Z_H_vals, y=Z_p_vals, fill=True, color="m", alpha=0.30, levels=25, ax=ax, label="Solution KDE")

    # Canonical markers
    for name, ZH, Zp, m, c in canon_pts:
        ax.plot(ZH, Zp, marker=m, ms=9, linestyle='none', color=c, label=name.capitalize())

    # Extra sequence markers
    for name, ZH, Zp, m, c in extra_pts:
        ax.plot(ZH, Zp, marker=m, ms=10, linestyle='none', color=c, label=name)

    # Two real overlays
    for ZH, Zp, lab, col, mkr, outline in real_sets:
        if ZH.size == 0: continue
        if outline:
            ax.scatter(ZH, Zp, marker=mkr, s=36, facecolors="none", edgecolors=col,
                       linewidths=1.0, alpha=0.95, label=f"{lab} (n={ZH.size})", zorder=3)
        else:
            ax.scatter(ZH, Zp, marker=mkr, s=30, c=col, alpha=0.85,
                       edgecolors="black", linewidths=0.5, label=f"{lab} (n={ZH.size})", zorder=4)

    # Cosmetics
    ax.axvline(0, lw=1, alpha=0.5)
    ax.axhline(0, lw=1, alpha=0.5)
    ax.set_xlabel("Z(H)  (vs. random at fixed L,k)")
    ax.set_ylabel("Z(p)  (vs. random at fixed L,k)")
    ax.set_title("Arbitrary vs. Real Sequence Feasibility")
    ax.legend(loc="best", frameon=True, fontsize=10)

    plt.tight_layout()
    plt.savefig(OUTFILE, dpi=300)
    plt.show()

    progress(f"[ok] saved → {OUTFILE}")
