#!/bin/bash
#SBATCH --time 24:00:00
#SBATCH --nodes=6
#SBATCH -p batch
#SBATCH --job-name spaceFISTA_0
#SBATCH -c 32

module load matlab
index=$( printf '%i' ${SLURM_ARRAY_TASK_ID} )
matlab -nodisplay -nodesktop -nosplash -r "ringFit2D_cyclic_main_0; exit"
#matlab -nodisplay -nodesktop -nosplash -r "addpath('../CHESS-Research/MATLAB'); ringFit2D_cyclic_main_0; exit"

