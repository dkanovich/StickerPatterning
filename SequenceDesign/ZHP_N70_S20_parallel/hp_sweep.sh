#!/bin/bash
#SBATCH --job-name=feas_ZH_ZP
#SBATCH --array=0-97%100          # 7×14=98 jobs; run up to 100 at once
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8          # <= set how many restarts to run in parallel
#SBATCH --mem=1G
#SBATCH --time=3-00:00:00
#SBATCH --partition=shakhnovich
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

module load Anaconda2/2019.10-fasrc01            
conda activate activate H_pIdx                      

mkdir -p logs results

# Avoid BLAS oversubscription inside workers
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# ── grid specification (Z-space) ────────────────────────────────────────
H_MIN=-1.0
H_MAX=5.0
H_N=7                      # number of H points (inclusive)

P_MIN=-4.0
P_MAX=9.0
P_N=14                      # number of p points (inclusive)

# Step sizes (floating point bc math)
H_STEP=$(echo "scale=12; ($H_MAX - $H_MIN) / ($H_N - 1)" | bc -l)
P_STEP=$(echo "scale=12; ($P_MAX - $P_MIN) / ($P_N - 1)" | bc -l)

# Map array index → (H_IDX, P_IDX)
GRID_WIDTH=$P_N
H_IDX=$(( SLURM_ARRAY_TASK_ID / GRID_WIDTH ))
P_IDX=$(( SLURM_ARRAY_TASK_ID % GRID_WIDTH ))

# Guard in case someone changes --array range
if (( H_IDX < 0 || H_IDX >= H_N || P_IDX < 0 || P_IDX >= P_N )); then
  echo "Array index $SLURM_ARRAY_TASK_ID out of range (H_IDX=$H_IDX, P_IDX=$P_IDX)."
  exit 1
fi

# Targets (keep 3 decimals for readability)
H_TARGET=$(printf "%.3f" "$(echo "$H_MIN + $H_IDX * $H_STEP" | bc -l)")
P_TARGET=$(printf "%.3f" "$(echo "$P_MIN + $P_IDX * $P_STEP" | bc -l)")
echo "Task $SLURM_ARRAY_TASK_ID ⇒ Z(H)=$H_TARGET  Z(p)=$P_TARGET"

# ── optimizer settings ─────────────────────────────────────────────────
NSTICKER=20
LENGTH=70
SEED=$SLURM_ARRAY_TASK_ID

# Parallel restarts inside the task
WORKERS=${SLURM_CPUS_PER_TASK:-1}
RESTARTS=${RESTARTS:-$WORKERS}     # by default, run one restart per core

OUTDIR=$PWD/results
OUTFILE=$OUTDIR/H${H_TARGET}_p${P_TARGET}.json

# ── run ────────────────────────────────────────────────────────────────
srun --cpu-bind=cores python3 seq_hp_opt.py \
  --nsticker "$NSTICKER" --length "$LENGTH" \
  --htarget "$H_TARGET" --ptarget "$P_TARGET" \
  --seed "$SEED" \
  --workers "$WORKERS" --restarts "$RESTARTS" \
  --out "$OUTFILE"
