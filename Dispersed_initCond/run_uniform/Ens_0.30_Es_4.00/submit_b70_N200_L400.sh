#!/bin/bash
#SBATCH --job-name=rv_Ens_0.30_Es_4.00
#SBATCH --array=1-5
#SBATCH -n  28
#SBATCH -t 3-00:00:00 # DD-HH:MM:SS
#SBATCH -p sapphire
#SBATCH --mem-per-cpu=400
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=davidkanovich@gmail.com

module load gcc/12.2.0-fasrc01 openmpi/4.1.4-fasrc01


# Auto-select input for this array task
INPUT_FILE="b70_N200_L400_run${SLURM_ARRAY_TASK_ID}.in"

srun -n $SLURM_NTASKS --mpi=pmix /n/home00/dkanovich/lammps-5Jun19/src/lmp_mpi -in $INPUT_FILE 
