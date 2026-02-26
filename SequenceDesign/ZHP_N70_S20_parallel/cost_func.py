# ──────────────────────────────────────────────────────────────────────
# Cost_function from JSONs (any H/p ranges; smooth, log-scaled coloring)
# ──────────────────────────────────────────────────────────────────────
import os, glob, json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, LogFormatterMathtext
from mpl_toolkits.axes_grid1 import make_axes_locatable

# === User settings ===
RESULT_DIR    = "results"
JSON_GLOB     = "*.json"
TARGET_KEYS_H = ["H_target", "h_target", "H_tgt", "h_tgt"]
TARGET_KEYS_P = ["p_target", "pIdx_target", "p_tgt", "pidx_target"]
COST_KEYS     = ["cost", "best_cost", "final_cost", "objective", "loss"]
NESTED_BUCKET = ["metrics", "result", "summary"]  # places cost might live

def first_present(d, keys):
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return None

def find_cost(d):
    v = first_present(d, COST_KEYS)
    if v is not None: return v
    for b in NESTED_BUCKET:
        if isinstance(d.get(b), dict):
            v = first_present(d[b], COST_KEYS)
            if v is not None: return v
    hist = d.get("history")
    if isinstance(hist, list) and hist:
        for rec in reversed(hist):
            if isinstance(rec, dict):
                v = first_present(rec, COST_KEYS)
                if v is not None: return v
    return None

def to_float(x):
    try:    return float(x)
    except: return None

def padded_limits(vals, frac=0.02, min_pad=1e-8):
    a = np.asarray(vals, float)
    a = a[np.isfinite(a)]
    lo, hi = float(a.min()), float(a.max())
    if lo == hi:
        span = max(1.0, abs(lo))
        lo, hi = lo - 0.5*span, hi + 0.5*span
    span = hi - lo
    pad = max(min_pad, frac*span)
    return lo - pad, hi + pad

# ---- load data → x,y,z ------------------------------------------------
paths = sorted(glob.glob(os.path.join(RESULT_DIR, JSON_GLOB)))
if not paths:
    raise FileNotFoundError(f"No JSON files found in {os.path.join(RESULT_DIR, JSON_GLOB)}")

xs, ys, zs = [], [], []
bad = 0
for pth in paths:
    try:
        with open(pth) as fh:
            d = json.load(fh)
    except Exception:
        bad += 1
        continue

    Ht = to_float(first_present(d, TARGET_KEYS_H))
    pt = to_float(first_present(d, TARGET_KEYS_P))
    c  = to_float(find_cost(d))

    if (Ht is not None and pt is not None and c is not None
        and np.isfinite(Ht) and np.isfinite(pt) and np.isfinite(c) and c > 0):
        xs.append(Ht); ys.append(pt); zs.append(c)
    else:
        bad += 1

if not xs:
    raise RuntimeError("No valid (H_target, p_target, cost>0) triples. Check your JSON keys.")

x = np.asarray(xs, float)
y = np.asarray(ys, float)
z = np.asarray(zs, float)

print(f"[info] loaded: {len(x)} points; skipped: {bad}")
print(f"[info] H range [{x.min():.3g}, {x.max():.3g}]  p range [{y.min():.3g}, {y.max():.3g}]  "
      f"cost range [{z.min():.3g}, {z.max():.3g}]")

# ---- plotting ----------------------------------------------------------
x_min, x_max = padded_limits(x)
y_min, y_max = padded_limits(y)

fig, ax = plt.subplots(figsize=(6, 6))
cm = matplotlib.colormaps.get_cmap('viridis')

# shared log norm (robust to outliers)
vmin = np.nanpercentile(z, 1)
vmax = np.nanpercentile(z, 99)
if not np.isfinite(vmin) or vmin <= 0: vmin = np.nanmin(z[z > 0])
if not np.isfinite(vmax):              vmax = np.nanmax(z)
if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
    vmin, vmax = z.min(), z.max()
    if vmin <= 0: vmin = np.nextafter(0, 1)  # tiny positive
norm = LogNorm(vmin=vmin, vmax=vmax)

# smooth field over scattered points
if len(x) >= 3:
    tri = mtri.Triangulation(x, y)
    ax.tricontourf(tri, z, levels=64, cmap=cm, norm=norm, alpha=0.85, zorder=1)
else:
    print("[warn] Not enough points for tricontourf; showing scatter only.")

# scatter (same norm so colors match)
scat = ax.scatter(x, y, c=z, cmap=cm, norm=norm, s=14, edgecolors='none', alpha=0.8, zorder=2)

# colorbar with decade ticks
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)
cbar = fig.colorbar(scat, cax=cax, orientation='vertical')
cbar.locator = LogLocator(base=10, numticks=7)
cbar.formatter = LogFormatterMathtext(base=10)
cbar.update_ticks()

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_xlabel('H')
ax.set_ylabel('pIdx')
ax.set_title('Cost_function (targets colored by cost, log scale)')
plt.tight_layout()
plt.savefig('cost_func_from_json.svg', dpi=300)
plt.show()
