#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NormMetric.py (Z-only)
----------------------
Provides Z-score normalization for dispersion metrics:
  • normalised_H(seq, sticker='S', lookup_path='H_Z_lookup.pkl')  → Z(H)
  • normalised_P(seq, sticker='S', lookup_path='P_Z_lookup.pkl')  → Z(p)
Both functions robustly handle lookup tables keyed by L or (L,k).
"""
import pickle, pathlib
from functools import lru_cache
from typing import Any, Dict, Tuple
from PatternMetric import StickerPattern

@lru_cache(maxsize=None)
def _load_lookup(path: str | pathlib.Path) -> Dict[Any, Any]:
    p = pathlib.Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Lookup file not found: {p}")
    with p.open("rb") as fh:
        return pickle.load(fh)

def _get_mu_sigma(entry_L: Any, k: int) -> Tuple[float, float]:
    """Extract (mu, sigma) from either (mu, sigma) or {k: (mu,sigma)}."""
    # length-only form
    if isinstance(entry_L, (tuple, list)):
        if len(entry_L) < 2:
            raise ValueError(f"Bad lookup tuple: {entry_L}")
        mu, sig = float(entry_L[0]), float(entry_L[1])
        return mu, max(sig, 1e-12)
    # length+count form
    if isinstance(entry_L, dict):
        # try int, then str, then numeric match
        if k in entry_L:
            mu, sig = entry_L[k]
        elif str(k) in entry_L:
            mu, sig = entry_L[str(k)]
        else:
            mu = sig = None
            for kk, val in entry_L.items():
                try:
                    if int(kk) == int(k):
                        mu, sig = val
                        break
                except Exception:
                    continue
            if mu is None:
                raise KeyError(f"Sticker count k={k} not found in length+count table.")
        return float(mu), max(float(sig), 1e-12)
    raise TypeError(f"Unexpected LOOKUP[L] type: {type(entry_L)}")

def _lookup_stats(LOOKUP: Dict[Any, Any], L: int, k: int) -> Tuple[float, float]:
    """Find (mu, sigma) for a given length L and sticker count k."""
    # direct int key
    entry_L = LOOKUP.get(L)
    # direct str key
    if entry_L is None:
        entry_L = LOOKUP.get(str(L))
    # numeric match across weird keys
    if entry_L is None:
        for LL, val in LOOKUP.items():
            try:
                if int(LL) == int(L):
                    entry_L = val
                    break
            except Exception:
                continue
    if entry_L is None:
        raise KeyError(f"Length L={L} not found in lookup.")
    return _get_mu_sigma(entry_L, k)

def normalised_H(seq: str, sticker: str = 'S', lookup_path: str = 'H_Z_lookup_100k.pkl') -> float:
    """Return Z-score for entropy H (frameshift-averaged over bL in PatternMetric)."""
    H_raw = StickerPattern(seq, sticker).calc_H_ave()
    L = len(seq)
    k = seq.count(sticker)
    LOOKUP = _load_lookup(lookup_path)
    mu, sig = _lookup_stats(LOOKUP, L, k)
    return 0.0 if sig <= 0 else (H_raw - mu) / sig

def normalised_P(seq: str, sticker: str = 'S', lookup_path: str = 'P_Z_lookup_100k.pkl') -> float:
    """Return Z-score for patchiness p (adjacent-sticker fraction)."""
    p_raw = StickerPattern(seq, sticker).process_patch()
    L = len(seq)
    k = seq.count(sticker)
    LOOKUP = _load_lookup(lookup_path)
    mu, sig = _lookup_stats(LOOKUP, L, k)
    return 0.0 if sig <= 0 else (p_raw - mu) / sig
