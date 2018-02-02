#!/bin/bash
#SBATCH --time 24:00:00
#SBATCH --nodes=2
#SBATCH --job-name spaceFISTAtest

module load matlab
index=$( printf '%i' ${SLURM_ARRAY_TASK_ID} )
matlab -nodisplay -nodesktop -nosplash -r "addpath('../CHESS-Research/MATLAB'); ringFit2D_cyclic_main; exit"

