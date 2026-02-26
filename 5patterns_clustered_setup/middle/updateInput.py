import re, sys  
import numpy as np  

def writeInputFile(file, L=500):
    lines = []
    with open(file, 'r') as tf:
        lines = tf.readlines()
        for i, line in enumerate(lines):
            if re.search('variable fName string', line):
                lines[i] = f'variable fName string b70_N200_L{L} \n'

            if re.search('read_data', line):
                lines[i] = f'read_data b70_N200_L{L}.data extra/special/per/atom 10 \n'

            if re.search('fix CV_Rg all colvars', line):
                lines[i] = f'fix CV_Rg all colvars N200_Rg_L{L}.colvars' + ' output ${fName} \n'

    
    with open(f'b70_N200_L{L}' + '.in', 'w') as tmf:
        tmf.writelines(lines)

def writeSubmitScript(L):
    infile = f'b70_N200_L{L}' + '.in'

    with open(f'submit_b70_N200_L{L}.sh', 'w') as tf:
        tf.write('#!/bin/bash\n') 
        tf.write(f'#SBATCH --job-name=L{L}_b70\n')
        tf.write(f'#SBATCH -n  14\n')
        tf.write('#SBATCH -t 3-00:00:00 # DD-HH:MM:SS\n')
        tf.write('#SBATCH -p shared\n')
        tf.write(f'#SBATCH --mem-per-cpu=400\n')
        tf.write('#SBATCH -o %x_%j.out\n')
        tf.write('#SBATCH -e %x_%j.err\n')
        tf.write('#SBATCH --mail-type=END\n')
        tf.write('#SBATCH --mail-user=davidkanovich@gmail.com\n\n')
        tf.write('module load gcc/12.2.0-fasrc01 openmpi/4.1.4-fasrc01\n\n')
        tf.write(f'srun -n $SLURM_NTASKS --mpi=pmix /n/home00/dkanovich/lammps-5Jun19/src/lmp_mpi -in {infile} \n')
        
        

_, file, L = sys.argv 
writeInputFile(file, L) 
writeSubmitScript(L)

