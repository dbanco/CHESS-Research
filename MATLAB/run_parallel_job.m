function run_parallel_job(num_jobs)
parfor idx = 1:num_jobs
    slurm_batch_wrapper(idx);
end