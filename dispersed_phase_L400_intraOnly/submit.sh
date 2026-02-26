#!/bin/bash
###############################################################################
# submit_all_runs.sh – launch many LAMMPS jobs in nested run_<pattern>/COND dirs
###############################################################################
#SBATCH --job-name=launch_L400
#SBATCH --time=00:10:00          # enough for dozens of sbatch calls
#SBATCH --partition=sapphire
#SBATCH --output=launch_%j.out
#SBATCH --error=launch_%j.err
# (No CPUs or memory requested – the work here is just submitting jobs)

# ---------------------------------------------------------------------------
# 1) CONFIGURATION  –––>  >>> EDIT THESE <<<  –––––––––––––––––––––------------
# ---------------------------------------------------------------------------
BASE_DIR=~/2024/dispersed_phase_L400_intraOnly          # top-level folder

PATTERNS=(term middle extreme patch uniform)       # run_<pattern>
CONDITIONS=(Ens_0.30_Es_4.00 Ens_0.30_Es_6.00 Ens_0.30_Es_8.00)

SUBMIT_SCRIPT=submit_b70_N200_L400.sh              # script to sbatch
# ---------------------------------------------------------------------------

echo "Launching jobs from $(hostname) on $(date)"
echo "Base dir: $BASE_DIR"
echo

for pat in "${PATTERNS[@]}"; do
    RUN_DIR="${BASE_DIR}/run_${pat}"
    if [[ ! -d $RUN_DIR ]]; then
        echo "[WARN] $RUN_DIR absent – skipping pattern '${pat}'"
        continue
    fi

    for cond in "${CONDITIONS[@]}"; do
        JOB_DIR="${RUN_DIR}/${cond}"
        if [[ ! -d $JOB_DIR ]]; then
            echo "[WARN] $JOB_DIR absent – skipping condition '${cond}'"
            continue
        fi

        echo "→  Submitting inside  ${JOB_DIR}"
        (
            cd "$JOB_DIR" || { echo "   cd failed – skipped"; exit; }
            if [[ -f $SUBMIT_SCRIPT ]]; then
                sbatch "$SUBMIT_SCRIPT"
            else
                echo "   [ERR] $SUBMIT_SCRIPT not found – skipped"
            fi
        )
    done
done

echo
echo "Done! All available jobs have been queued."
