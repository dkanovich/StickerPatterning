"""
Map (H_target, p_target) → (H_seq, p_seq) for all optimisation jobs
===================================================================

• Expects a folder called  results/   in the current working directory,
  containing JSON files written by seq_hp_opt.py.  Each JSON must have
  the keys:  H_target, p_target, H_seq, p_seq.        (seed, best_seq,
  final_cost are ignored here.)

• Produces the figure displayed on-screen; uncomment the savefig() line
  if you want to save it.
"""

import json, glob, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ────────────────────────────── I/O ───────────────────────────────────
RESULT_DIR = "results"                       # change if you renamed it
json_paths = glob.glob(os.path.join(RESULT_DIR, "*.json"))
if not json_paths:
    raise FileNotFoundError(f"No JSON files found in {RESULT_DIR}/")

def first_present(d, keys):
    """Return value for the first key that exists in dict d (allow 0/False)."""
    for k in keys:
        if k in d:
            return d[k]
    return None

records = []
for path in json_paths:
    with open(path) as fh:
        data = json.load(fh)

    # Robust to naming + preserves zeros
    H_tgt = first_present(data, ["H_target", "h_target"])
    p_tgt = first_present(data, ["p_target", "pIdx_target", "p_tgt", "pidx_target"])
    H_sol = first_present(data, ["H_seq", "h_seq", "H"])
    p_sol = first_present(data, ["p_seq", "p"])

    # Optional: cast to float if present; keep None if missing
    try:    H_tgt = float(H_tgt) if H_tgt is not None else None
    except: pass
    try:    p_tgt = float(p_tgt) if p_tgt is not None else None
    except: pass
    try:    H_sol = float(H_sol) if H_sol is not None else None
    except: pass
    try:    p_sol = float(p_sol) if p_sol is not None else None
    except: pass

    records.append((H_tgt, p_tgt, H_sol, p_sol))

df = pd.DataFrame(records, columns=["H_target", "p_target", "H", "p"])\
       .sort_values(["H_target", "p_target"]).reset_index(drop=True)

# Drop rows with any missing required values (avoid plotting NaNs)
df = df.dropna(subset=["H_target", "p_target", "H", "p"])

# ──────────────────────────── PLOT ────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 6))

# Density of *solutions* (pink contour)
if len(df) > 1:
    sns.kdeplot(x=df["H"], y=df["p"], fill=True, color="m", alpha=0.3,
                label="Solution", ax=ax)

# Scatter raw points
ax.scatter(df["H_target"], df["p_target"], color="c", s=20, label="Target", alpha=0.6)
ax.scatter(df["H"],        df["p"],        color="m", s=20, label="Solution", alpha=0.6)

# Quiver arrows
U = df["H"] - df["H_target"]
V = df["p"] - df["p_target"]
ax.quiver(df["H_target"], df["p_target"], U, V,
          angles="xy", scale_units="xy", scale=1,
          color="black", width=0.002, alpha=0.8)

# Cosmetics
ax.set_xlabel("H")
ax.set_ylabel("p")
ax.set_title("Target to solution mapping")
ax.legend(loc="lower left", frameon=True)

# Save if desired
plt.savefig("./mapping_from_json2.svg")
plt.tight_layout()
plt.show()
