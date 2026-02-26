#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Most-common solution (Z(H), Z(P)) pairs across optimization runs
with exemplar sequence (lowest loss) per solution bin.
"""

import os, glob, json, math
from collections import Counter, defaultdict
from pathlib import Path

# ========================= CONFIG =========================
RESULTS_DIRS = [
    "results",          # add more folders if needed
]
DECIMALS = 2           # binning precision for Z-values
INCLUDE_TARGET_LIST = False  # keep False unless you want all target bins per solution
OUTPUT_TXT = "most_common_Z_solutions.txt"
PRINT_SKIP_STATS = True
# ==========================================================

# ---- Helpers -----------------------------------------------------------
def first_present(d, keys):
    for k in keys:
        if k in d:
            return d[k]
    return None

def to_float_or_none(x):
    try:
        return float(x)
    except Exception:
        return None

def invnorm_from_P(P):
    try:
        if P is None:
            return None
        P = float(P)
        if not (0.0 < P < 1.0):
            return None
        return math.sqrt(2.0) * math.erfinv(2.0*P - 1.0)
    except Exception:
        return None

def extract_Z_pair(d, zh_keys, zp_keys, ph_keys=None, pp_keys=None):
    zH = to_float_or_none(first_present(d, zh_keys))
    zP = to_float_or_none(first_present(d, zp_keys))
    if zH is None and ph_keys:
        zH = invnorm_from_P(to_float_or_none(first_present(d, ph_keys)))
    if zP is None and pp_keys:
        zP = invnorm_from_P(to_float_or_none(first_present(d, pp_keys)))
    return zH, zP

def round_pair(zH, zP, decimals):
    if zH is None or zP is None:
        return None
    if any(map(math.isnan, [zH, zP])):
        return None
    return (round(zH, decimals), round(zP, decimals))

def find_json_paths(dirs):
    paths = []
    for root in dirs:
        root_path = Path(root)
        if not root_path.exists():
            continue
        paths.extend([str(p) for p in root_path.rglob("*.json")])
    return sorted(paths)

# ======= KEY MAPS (accepts H_* / p_* and Z_* names) =====================
SOLUTION_ZH_KEYS = [
    "Z_H_seq","Z_H","ZH_seq","ZH","Z(H)_seq","Z(H)","Z_H_solution",
    "H_seq","H",
]
SOLUTION_ZP_KEYS = [
    "Z_p_seq","Z_p","ZP_seq","ZP","Z(p)_seq","Z(p)","Z_p_solution",
    "p_seq","p",
]
TARGET_ZH_KEYS = ["Z_H_target","ZH_target","Z(H)_target","H_target","h_target"]
TARGET_ZP_KEYS = ["Z_p_target","ZP_target","Z(p)_target","p_target","pIdx_target","p_tgt","pidx_target"]

SOLUTION_PH_KEYS = ["P_H_seq","P_H","PH_seq","PH","P(H)"]
SOLUTION_PP_KEYS = ["P_p_seq","P_p","PP_seq","PP","P(p)"]
TARGET_PH_KEYS   = ["P_H_target","PH_target","P(H)_target"]
TARGET_PP_KEYS   = ["P_p_target","PP_target","P(p)_target"]

# Optional fields
FINAL_COST_KEYS = ["final_cost","loss","cost","best_cost"]
BEST_SEQ_KEYS   = ["best_seq","sequence","seq"]

# ---- Main --------------------------------------------------------------
def main():
    json_paths = find_json_paths(RESULTS_DIRS)
    if not json_paths:
        raise FileNotFoundError(f"No JSON files found under: {RESULTS_DIRS}")

    counts = Counter()                 # solution bin -> count
    sol_to_targets = defaultdict(set)  # solution bin -> set of target bins (rounded)

    # For exemplar: solution bin -> dict with best (lowest) cost entry
    exemplars = {}  # (ZH_sol_bin, ZP_sol_bin) -> {"cost": float, "seq": str, "tgt_bin": (zHt,zPt)}

    total = used = 0
    skip_io = skip_no_sol = skip_nan = 0

    for p in json_paths:
        total += 1
        try:
            data = json.load(open(p, "r"))
        except Exception:
            skip_io += 1
            continue

        d = data["data"] if isinstance(data, dict) and isinstance(data.get("data"), dict) else data

        # Extract solution & target Z-pairs
        zH_sol, zP_sol = extract_Z_pair(
            d, SOLUTION_ZH_KEYS, SOLUTION_ZP_KEYS, SOLUTION_PH_KEYS, SOLUTION_PP_KEYS
        )
        if zH_sol is None or zP_sol is None:
            skip_no_sol += 1
            continue
        if any(map(lambda x: isinstance(x, float) and math.isnan(x), [zH_sol, zP_sol])):
            skip_nan += 1
            continue

        sol_bin = round_pair(zH_sol, zP_sol, DECIMALS)
        if sol_bin is None:
            skip_nan += 1
            continue

        counts[sol_bin] += 1
        used += 1

        zH_tgt, zP_tgt = extract_Z_pair(
            d, TARGET_ZH_KEYS, TARGET_ZP_KEYS, TARGET_PH_KEYS, TARGET_PP_KEYS
        )
        tgt_bin = round_pair(zH_tgt, zP_tgt, DECIMALS)
        if tgt_bin is not None:
            sol_to_targets[sol_bin].add(tgt_bin)

        # Exemplar tracking (lowest cost)
        cost = to_float_or_none(first_present(d, FINAL_COST_KEYS))
        seq  = first_present(d, BEST_SEQ_KEYS)
        # if missing cost, treat as +inf so any real cost will win
        cost_val = cost if cost is not None else float("inf")

        ex = exemplars.get(sol_bin)
        if ex is None or cost_val < ex["cost"]:
            exemplars[sol_bin] = {
                "cost": cost_val,
                "seq": seq if isinstance(seq, str) else "<no_sequence_in_json>",
                "tgt_bin": tgt_bin,   # may be None if target missing
            }

    if used == 0:
        msg = (
            "Parsed JSONs, but no valid (Z_H, Z_P) solutions found.\n"
            f"Total files = {total}\n"
            f"  I/O/JSON errors        : {skip_io}\n"
            f"  Missing solution fields: {skip_no_sol}\n"
            f"  NaN/invalid values     : {skip_nan}\n"
            "Check that your files contain H_seq/p_seq or Z_* fields."
        )
        raise RuntimeError(msg)

    # Write ranked report
    rows = counts.most_common()
    with open(OUTPUT_TXT, "w") as out:
        out.write("Most-common solution (Z(H), Z(P)) pairs with exemplar sequences\n")
        out.write("================================================================\n")
        out.write(f"Searched folders : {', '.join(RESULTS_DIRS)}\n")
        out.write(f"Total JSON files : {total}\n")
        out.write(f"Valid solutions  : {used}\n")
        out.write(f"Skipped (I/O)    : {skip_io}\n")
        out.write(f"Skipped (missing): {skip_no_sol}\n")
        out.write(f"Skipped (NaN)    : {skip_nan}\n")
        out.write(f"Binning (decimals): {DECIMALS}\n\n")

        out.write(f"{'Rank':>4}  {'Count':>7}  {'Share':>7}  {'ZH_sol':>8}  {'ZP_sol':>8}  {'best_cost':>12}\n")
        out.write("-"*70 + "\n")

        for i, ((zh, zp), cnt) in enumerate(rows, start=1):
            share = 100.0 * cnt / used
            ex = exemplars.get((zh, zp), {"cost": float("inf"), "seq": "<none>", "tgt_bin": None})
            best_cost = ex["cost"]
            out.write(f"{i:>4}  {cnt:>7}  {share:6.1f}%  {zh:8.2f}  {zp:8.2f}  {best_cost:12.6g}\n")

            # Exemplar details (sequence + the target it came from)
            tgt_bin = ex.get("tgt_bin")
            if tgt_bin is not None:
                tzh, tzp = tgt_bin
                out.write(f"      exemplar_from_target: ZH_target={tzh:.{DECIMALS}f}, ZP_target={tzp:.{DECIMALS}f}\n")
            else:
                out.write(f"      exemplar_from_target: (missing in JSON)\n")

            out.write("      exemplar_sequence   : ")
            out.write(ex["seq"] if isinstance(ex["seq"], str) else "<invalid sequence type>")
            out.write("\n")

            if INCLUDE_TARGET_LIST:
                tgts = sorted(sol_to_targets.get((zh, zp), []))
                if tgts:
                    out.write("      all_target_bins     : ")
                    out.write(", ".join(f"({tzh:.2f},{tzp:.2f})" for tzh, tzp in tgts))
                    out.write("\n")
            out.write("\n")

        out.write("-"*70 + "\n")
        out.write(f"Total valid = {used} (100.0%)\n")

    if PRINT_SKIP_STATS:
        print(f"[info] files: total={total}, used={used}, skip_io={skip_io}, "
              f"skip_missing={skip_no_sol}, skip_nan={skip_nan}")
    print(f"[done] Wrote → {OUTPUT_TXT}")
    print(f"[info] Unique solution bins: {len(rows)} (DECIMALS={DECIMALS})")

if __name__ == "__main__":
    main()
