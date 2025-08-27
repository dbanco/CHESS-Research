function error_stats = compute_error_stats(objectives, sig_ind)
% COMPUTE_ERROR_STATS Computes mean and std of error metrics for given indices.
%   Returns a struct containing all output values.
%
%   INPUTS:
%       objectives_indep - cell array of structs with fields: error, rel_error, true_error
%       sig_ind - array of indices to evaluate
%
%   OUTPUT:
%       error_stats - struct with fields:
%           avg_error, std_error,
%           avg_rel_error, std_rel_error,
%           avg_true_error, std_true_error,
%           avg_Derror, std_Derror

    % Preallocate full-size vectors
    max_ind = numel(sig_ind);
    error_stats.avg_error      = zeros(max_ind, 1);
    error_stats.std_error      = zeros(max_ind, 1);
    error_stats.avg_rel_error  = zeros(max_ind, 1);
    error_stats.std_rel_error  = zeros(max_ind, 1);
    error_stats.avg_true_error = zeros(max_ind, 1);
    error_stats.std_true_error = zeros(max_ind, 1);
    error_stats.avg_Derror     = zeros(max_ind, 1);
    error_stats.std_Derror     = zeros(max_ind, 1);
    error_stats.avg_of_penalty = zeros(max_ind, 1);
    error_stats.std_of_penalty = zeros(max_ind, 1);

    % Compute statistics
    i = 1;
    for n = sig_ind
        error_stats.avg_error(i)      = mean(objectives{i}.error,'omitnan');
        error_stats.std_error(i)      = std(objectives{i}.error,'omitnan');
        error_stats.avg_rel_error(i)  = mean(objectives{i}.rel_error,'omitnan');
        error_stats.std_rel_error(i)  = std(objectives{i}.rel_error,'omitnan');
        error_stats.avg_true_error(i) = mean(objectives{i}.true_error,'omitnan');
        error_stats.std_true_error(i) = std(objectives{i}.true_error,'omitnan');
        error_stats.avg_Derror(i)     = mean(objectives{i}.D_error,'omitnan'); 
        error_stats.std_Derror(i)     = std(objectives{i}.D_error,'omitnan'); 
        error_stats.avg_of_penalty(i)     = mean(objectives{i}.of_penalty,'omitnan'); 
        error_stats.std_of_penalty(i)     = std(objectives{i}.of_penalty,'omitnan'); 
        error_stats.avg_hs_penalty(i)     = mean(objectives{i}.hs_penalty,'omitnan'); 
        error_stats.std_hs_penalty(i)     = std(objectives{i}.hs_penalty,'omitnan'); 
        error_stats.avg_ofhs_penalty(i)     = mean(objectives{i}.of_penalty+objectives{i}.hs_penalty,'omitnan'); 
        error_stats.std_ofhs_penalty(i)     = std(objectives{i}.of_penalty+objectives{i}.hs_penalty,'omitnan'); 
        error_stats.avg_log_penalty(i)     = mean(objectives{i}.log_penalty,'omitnan'); 
        error_stats.std_log_penalty(i)     = std(objectives{i}.log_penalty,'omitnan'); 
        i = i + 1;
    end
end
