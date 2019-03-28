% Initialize unregularized solution
cd '../../job_wass_parallel_init'
parfor idx = 1:n
    slurm_batch_wrapper(idx);
end

% 10 pairs of regularized iterations
for i = 1:10
    cd '../job_wass_parallel_ab'
    parfor idx = 1:n
        slurm_batch_wrapper(idx);
    end
    
    cd '../job_wass_parallel_ba'
    parfor idx = 1:n
        slurm_batch_wrapper(idx);
    end
end