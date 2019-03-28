function run_parallel_job(num_jobs)
parpool(32)
parfor idx = 1:num_jobs
    slurm_batch_wrapper(idx);
end
