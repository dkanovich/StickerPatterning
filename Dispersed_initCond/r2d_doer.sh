#!/usr/bin/env bash
# deploy_r2d.sh
# ─────────────
# Copy r2d_hpc.sh into every energy folder and launch it with sbatch.
# Adjust ROOT_DIR or the PATTERNS/ENERGIES arrays if your layout changes.

set -euo pipefail

# --- USER SETTINGS ----------------------------------------------------------
ROOT_DIR="."   # project root
MASTER_R2D="${ROOT_DIR}/r2d_hpc.sh"                  # master copy to copy from
PATTERNS=(middle term extreme uniform patch)
ENERGIES=(Ens_0.30_Es_4.00 Ens_0.30_Es_6.00 Ens_0.30_Es_8.00)
# ----------------------------------------------------------------------------

# Check the master script existss
if [[ ! -f "$MASTER_R2D" ]]; then
  echo "✖  Master r2d_hpc.sh not found at $MASTER_R2D" >&2
  exit 1
fi

echo "Deploying r2d_hpc.sh from $MASTER_R2D …"
jobs_launched=0

for patt in "${PATTERNS[@]}"; do
  for en in "${ENERGIES[@]}"; do
    energy_dir="${ROOT_DIR}/run_${patt}/${en}"
    if [[ ! -d "$energy_dir" ]]; then
      echo "  [skip] $energy_dir (missing)"
      continue
    fi

    # 1) copy the script
    cp -f "$MASTER_R2D" "$energy_dir/"        # -f to overwrite silently
    echo "  [copy] ${energy_dir#"$ROOT_DIR/"}"

    # 2) launch the job from inside the energy folder
    (
      cd "$energy_dir"
      sbatch r2d_hpc.sh
    )
    jobs_launched=$((jobs_launched + 1))
  done
done

echo "✔  Done — $jobs_launched sbatch jobs submitted."
