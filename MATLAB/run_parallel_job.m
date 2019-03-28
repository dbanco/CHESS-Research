% Run from appropriate job folder
addpath('../CHESS-Research/MATLAB')

parfor idx = 1:2
    slurm_batch_wrapper(idx);
end