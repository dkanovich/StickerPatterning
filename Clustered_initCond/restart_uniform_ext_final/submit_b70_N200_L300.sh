#!/bin/bash
#SBATCH --job-name=uniform_prod
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

# Get the input file for this task based on the array index
INPUT_FILE="b70_N200_L300_run${SLURM_ARRAY_TASK_ID}.in"

# Run the simulation using the specific input file
srun -n $SLURM_NTASKS --mpi=pmix /n/home00/dkanovich/lammps-5Jun19/src/lmp_mpi -in $INPUT_FILE