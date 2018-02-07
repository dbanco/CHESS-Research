#!/bin/bash
#SBATCH --time 24:00:00
#SBATCH --nodes=6
#SBATCH -p batch
#SBATCH --job-name spaceFISTA
#SBATCH -c 32
#SBATCH --mail-type=END
#SBATCH --mail-user=dbanco02@tufts.edu

module load matlab
#index=$( printf '%i' ${SLURM_ARRAY_TASK_ID} )
matlab -nodisplay -nodesktop -nosplash -r "ringFit2D_cyclic_main_test; exit"
#matlab -nodisplay -nodesktop -nosplash -r "addpath('../CHESS-Research/MATLAB'); ringFit2D_cyclic_main_0; exit"

