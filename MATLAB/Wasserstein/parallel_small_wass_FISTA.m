function parallel_small_wass_FISTA(n,jobDirName)
parpool(32)

dir1 = ['../job_',jobDirName,'_ab'];
dir2 = ['../job_',jobDirName,'_ba'];

% 10 pairs of regularized iterations
for i = 1:200
    
    cd(dir1)
    parfor idx = 1:n
        slurm_batch_wrapper(idx);
    end
    
    cd(dir2)
    parfor idx = 1:n
        slurm_batch_wrapper(idx);
    end
end
