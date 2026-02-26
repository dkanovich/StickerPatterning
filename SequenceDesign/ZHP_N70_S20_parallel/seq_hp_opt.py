#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
seq_hp_opt.py - optimize ONE (Z(H)_target, Z(p)_target) point.
Adds parallel multi-start via --restarts and --workers.
- Targets are interpreted as Z-scores (Z-only pipeline).
- Output keys remain backward-compatible: H_target/p_target (targets in Z),
  H_seq/p_seq (achieved Z), plus final_cost, cost_curve, best_seq, seed.
"""
import argparse
import json
import os
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict

import numpy as np

from SeqGen_Feasibility_space import SeqDesign
from NormMetric import normalised_H, normalised_P


def _one_run(seed: int, nsticker: int, length: int,
             zH_tgt: float, zP_tgt: float,
             nstep: int, steps: int, sticker: str) -> Dict:
    """Run one independent restart and return its result bundle."""
    print(f"[restart] seed={seed}  targets: Z(H)={zH_tgt:.3f}, Z(p)={zP_tgt:.3f}", flush=True)

    # Reproducibility per restart
    random.seed(seed)
    np.random.seed(seed)

    sd = SeqDesign(nsticker, length, zH_tgt, zP_tgt, sticker_char=sticker)
    best_seq, epochs, costs = sd.generate_seq(nStep=nstep, steps=steps)

    # Achieved Z metrics for the best sequence
    zH = float(normalised_H(best_seq, sticker=sticker))
    zP = float(normalised_P(best_seq, sticker=sticker))

    # Ensure final_cost corresponds to the reported best_seq
    final_cost = 0.5 * (zH - zH_tgt) ** 2 + 0.5 * (zP - zP_tgt) ** 2

    return {
        "seed": seed,
        "best_seq": best_seq,
        "zH": zH,
        "zP": zP,
        "final_cost": float(final_cost),
        "cost_curve": [float(c) for c in costs],
    }


def main():
    ap = argparse.ArgumentParser(description="Optimize a single Z(H), Z(p) target with parallel restarts.")
    ap.add_argument("--nsticker", type=int, required=True, help="Number of stickers in the sequence.")
    ap.add_argument("--length",   type=int, required=True, help="Total sequence length.")
    ap.add_argument("--htarget",  type=float, required=True, help="Target Z(H).")
    ap.add_argument("--ptarget",  type=float, required=True, help="Target Z(p).")
    ap.add_argument("--seed",     type=int, default=0, help="Base RNG seed.")
    ap.add_argument("--out",      type=str, required=True, help="Path to JSON output.")
    # New knobs for scaling and schedule
    ap.add_argument("--restarts", type=int, default=1, help="Independent restarts to run (best is kept).")
    ap.add_argument("--workers",  type=int, default=1, help="Parallel processes to use.")
    ap.add_argument("--nstep",    type=int, default=12, help="Annealing epochs per restart.")
    ap.add_argument("--steps",    type=int, default=20000, help="Swap attempts per epoch.")
    ap.add_argument("--sticker",  type=str, default="S", help="Sticker character used in the sequence.")
    args = ap.parse_args()

    # Avoid BLAS oversubscription in workers (harmless if libs absent)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    jobs = []
    base_seed = int(args.seed)
    restarts = max(1, int(args.restarts))
    workers  = max(1, int(args.workers))

    if workers > 1 and restarts > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = []
            for r in range(restarts):
                futs.append(ex.submit(
                    _one_run,
                    base_seed + r,
                    int(args.nsticker),
                    int(args.length),
                    float(args.htarget),
                    float(args.ptarget),
                    int(args.nstep),
                    int(args.steps),
                    str(args.sticker),
                ))
            for f in as_completed(futs):
                jobs.append(f.result())
    else:
        # Run sequentially (workers==1 or restarts==1)
        for r in range(restarts):
            jobs.append(_one_run(
                base_seed + r,
                int(args.nsticker),
                int(args.length),
                float(args.htarget),
                float(args.ptarget),
                int(args.nstep),
                int(args.steps),
                str(args.sticker),
            ))

    # Pick the best by final_cost
    best = min(jobs, key=lambda d: d["final_cost"])
    
    print(f"[best] seed={best['seed']}  Z(H)={best['zH']:.3f}  Z(p)={best['zP']:.3f}  "
      f"final_cost={best['final_cost']:.6g}", flush=True)

    out = {
        "seed":       best["seed"],
        "restarts":   restarts,
        "workers":    workers,
        "H_target":   float(args.htarget),   # Z(H) target
        "p_target":   float(args.ptarget),   # Z(p) target
        "H_seq":      best["zH"],            # achieved Z(H)
        "p_seq":      best["zP"],            # achieved Z(p)
        "final_cost": best["final_cost"],
        "cost_curve": best["cost_curve"],
        "best_seq":   best["best_seq"],
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
