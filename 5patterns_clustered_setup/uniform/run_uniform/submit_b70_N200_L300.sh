#!/bin/bash
#SBATCH --job-name=uniform
#SBATCH -n  28
#SBATCH -t 3-00:00:00 # DD-HH:MM:SS
#SBATCH -p sapphire
#SBATCH --mem-per-cpu=400
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=davidkanovich@gmail.com

module load gcc/12.2.0-fasrc01 openmpi/4.1.4-fasrc01

srun -n $SLURM_NTASKS --mpi=pmix /n/home00/dkanovich/lammps-5Jun19/src/lmp_mpi -in b70_N200_L300.in 
