function parallel_reg_awmv_FISTA(n)
parpool(32)
% Initialize unregularized solution
parfor idx = 1:n
   slurm_batch_wrapper(idx);
end

% 10 pairs of regularized iterations
for i = 1:200
    cd '../job_awmv_parallel_ab'
    parfor idx = 1:n
        slurm_batch_wrapper(idx);
    end
    
    cd '../job_awmv_parallel_ba'
    parfor idx = 1:n
        slurm_batch_wrapper(idx);
    end
end
