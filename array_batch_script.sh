#!/bin/bash
#SBATCH --time 24:00:00
#SBATCH --array=0-923

module load matlab
index=$( printf '%i' ${SLURM_ARRAY_TASK_ID} )
matlab -nodisplay -nodesktop -nosplash -r "addpath('../CHESS-Research/MATLAB'); slurm_batch_wrapper('$index'); exit"

