function slurm_write_bash( numJobs, jobDir, scriptFileName, indexString )
%slurm_write_bash Writes bash script to initiate batch job

scriptText = [...
'#!/bin/bash\n',...
'#SBATCH --time 01:00:00\n',...
'#SBATCH --mem-per-cpu=16000\n',...
'#SBATCH --array=',indexString,'\n',...
'#SBATCH -p batch\n',...
'#SBATCH --error=slurm_error.err\n',...
'#SBATCH --mail-type=END\n',...
'#SBATCH --mail-user=dbanco02@tufts.edu\n',...
'\n',...
'module load matlab/2020a\n',...
'index=$( printf ''%%i'' ${SLURM_ARRAY_TASK_ID} )\n',...
'matlab -nodisplay -nodesktop -nosplash -r',...
' "addpath(genpath(''../CHESS-Research/MATLAB''));',...
'slurm_batch_wrapper(''$index''); exit"\n',...
];

fileID = fopen(fullfile(jobDir,scriptFileName),'w');
fprintf(fileID,scriptText,numJobs);
fclose(fileID);

end
