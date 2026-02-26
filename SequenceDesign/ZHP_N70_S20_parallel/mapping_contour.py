#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Canonical + solution KDE in Z-space (X=Z(H), Y=Z(p))

• Reads results/*.json that already contain Z-scores:
      H_target, p_target, H_seq, p_seq
• Plots KDE of solutions exactly like your mapping script:
      sns.kdeplot(x=..., y=..., fill=True, color="m", alpha=0.3, label="Solution")
• Overlays canonical 5 patterns in Z-space (computed via lookups or MC fallback).
"""

import json, glob, os, pickle, pathlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ───────────────── I/O ─────────────────
RESULT_DIR    = "results"
json_paths    = glob.glob(os.path.join(RESULT_DIR, "*.json"))
if not json_paths:
    raise FileNotFoundError(f"No JSON files found in {RESULT_DIR}/")

# ───────────────── Load solutions (already Z) ─────────────────
def first_present(d, keys):
    for k in keys:
        if k in d:
            return d[k]
    return None

Z_H_vals, Z_p_vals = [], []
for path in json_paths:
    with open(path) as fh:
        data = json.load(fh)
    ZH = first_present(data, ["H_seq", "h_seq", "H"])  # already Z(H)
    Zp = first_present(data, ["p_seq", "p"])           # already Z(p)
    if ZH is None or Zp is None:
        continue
    try:
        Z_H_vals.append(float(ZH))
        Z_p_vals.append(float(Zp))
    except Exception:
        continue

Z_H_vals = np.asarray(Z_H_vals, dtype=float)
Z_p_vals = np.asarray(Z_p_vals, dtype=float)
n_used = len(Z_H_vals)

# ───────────────── Canonical patterns (compute Z for display) ─────────────────
# Only needed for the five markers. If lookups absent, MC fallback is used.
LOOKUP_PATH_H = pathlib.Path("H_Z_lookup_100k.pkl")
LOOKUP_PATH_P = pathlib.Path("P_Z_lookup_100k.pkl")
LOOKUP_H = pickle.loads(LOOKUP_PATH_H.read_bytes()) if LOOKUP_PATH_H.exists() else None
LOOKUP_P = pickle.loads(LOOKUP_PATH_P.read_bytes()) if LOOKUP_PATH_P.exists() else None

from numpy import array, log, mean
import numpy.ma as ma

def _calc_H_from_counts(counts):
    p_arr = array(counts) / sum(counts)
    return float(np.sum(-1 * p_arr * ma.log(p_arr).filled(0)))

def _H_once(seq: str, aa: str, bL: int, offset: int) -> float:
    N = len(seq)
    counts = []
    for j in range(0, N, bL):
        block = seq[j + offset : j + offset + bL]
        if len(block) < bL:
            block += seq[: bL - len(block)]
        counts.append(block.count(aa))
    return _calc_H_from_counts(counts) / log(len(counts))

def H_raw(seq: str, aa: str) -> float:
    return float(mean([mean([_H_once(seq, aa, bL, off) for off in range(bL)])
                       for bL in (4, 5)]))

def patchiness(seq: str, aa: str) -> float:
    N = seq.count(aa)
    if N <= 1: return 0.0
    adj = sum(1 for i in range(len(seq)-1) if seq[i]==aa and seq[i+1]==aa)
    return adj/(N-1)

def _lookup_mu_sigma(lookup, L, k):
    if lookup is None: return None
    entry_L = lookup.get(L, lookup.get(str(L)))
    if not isinstance(entry_L, dict): return None
    entry_k = entry_L.get(k, entry_L.get(str(k)))
    if entry_k is None: return None
    mu, sig = float(entry_k[0]), float(entry_k[1])
    return mu, max(sig, 1e-12)

# MC fallback for canonical only
def _mc_baseline(L, k, aa, n_ref=5000, seed=123):
    rng = np.random.default_rng(seed)
    base = np.array([aa]*k + ['2']*(L-k), dtype='<U1')
    H_vals = np.empty(n_ref, dtype=np.float32)
    p_vals = np.empty(n_ref, dtype=np.float32)
    for i in range(n_ref):
        rng.shuffle(base)
        s = ''.join(base)
        H_vals[i] = H_raw(s, aa)
        p_vals[i] = patchiness(s, aa)
    mu_H, sig_H = float(H_vals.mean()), float(H_vals.std(ddof=0))
    mu_p, sig_p = float(p_vals.mean()), float(p_vals.std(ddof=0))
    return mu_H, max(sig_H,1e-12), mu_p, max(sig_p,1e-12)

def to_Z_from_raw(H_val, p_val, L, k, aa='1'):
    mu_sig_H = _lookup_mu_sigma(LOOKUP_H, L, k)
    mu_sig_P = _lookup_mu_sigma(LOOKUP_P, L, k)
    if mu_sig_H is None or mu_sig_P is None:
        mu_H, sH, mu_p, sP = _mc_baseline(L, k, aa)
    else:
        mu_H, sH = mu_sig_H
        mu_p, sP = mu_sig_P
    return (H_val - mu_H)/sH, (p_val - mu_p)/sP

sequences = {
    "uniform": "2212212221221222122122212212221221222122122212212221221222122122212212",
    "term":    "1111111111111111111122222222222222222222222222222222222222222222222222",
    "middle":  "2222222222222222222222222111111111111111111112222222222222222222222222",
    "extreme": "1111111111222222222222222222222222222222222222222222222222221111111111",
    "patch":   "2222211112222222222111122222222221111222222222211112222222222111122222",
}
COLORS  = dict(term="green", uniform="red", middle="blue", extreme="purple", patch="orange")
MARKERS = dict(uniform='o', term='s', middle='^', extreme='D', patch='v')

canon_pts = []
for name, seq in sequences.items():
    L = len(seq)
    k = seq.count('1')
    H = H_raw(seq, '1')
    p = patchiness(seq, '1')
    ZH, Zp = to_Z_from_raw(H, p, L, k, '1')
    canon_pts.append((name, ZH, Zp))

# ───────────────── Plot ─────────────────
fig, ax = plt.subplots(figsize=(6, 6))

# KDE EXACTLY like your script
if n_used > 1:
    sns.kdeplot(x=Z_H_vals, y=Z_p_vals, fill=True, color="m", alpha=0.3,
                label="Solution", ax=ax)

# Optional: scatter the same solution points
if n_used > 0:
    ax.scatter(Z_H_vals, Z_p_vals, color="m", s=20, alpha=0.6, label="Solution pts")

# Canonical markers
for name, ZH, Zp in canon_pts:
    ax.plot(ZH, Zp, marker=MARKERS[name], ms=9, linestyle='none',
            color=COLORS[name], label=name.capitalize())

# Crosshairs & cosmetics
ax.axvline(0, lw=1, alpha=0.5)
ax.axhline(0, lw=1, alpha=0.5)
ax.set_xlabel("Z(H)  (vs. random at fixed L,k)")
ax.set_ylabel("Z(p)  (vs. random at fixed L,k)")
ax.set_title("Canonical + solutions in Z-space (X=Z(H), Y=Z(p))")
ax.legend(loc="best", frameon=True)

plt.tight_layout()
plt.savefig("canonical_ZH_vs_Zp_KDE_from_json.svg")
plt.show()

print(f"[info] solutions plotted: n={n_used}")
