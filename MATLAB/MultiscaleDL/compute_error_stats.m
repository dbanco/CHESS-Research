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
    max_ind = max(sig_ind);
    error_stats.avg_error      = zeros(max_ind, 1);
    error_stats.std_error      = zeros(max_ind, 1);
    error_stats.avg_rel_error  = zeros(max_ind, 1);
    error_stats.std_rel_error  = zeros(max_ind, 1);
    error_stats.avg_true_error = zeros(max_ind, 1);
    error_stats.std_true_error = zeros(max_ind, 1);
    error_stats.avg_Derror     = zeros(max_ind, 1);
    error_stats.std_Derror     = zeros(max_ind, 1);

    % Compute statistics
    for n = sig_ind
        error_stats.avg_error(n)      = mean(objectives{n}.error);
        error_stats.std_error(n)      = std(objectives{n}.error);
        error_stats.avg_rel_error(n)  = mean(objectives{n}.rel_error);
        error_stats.std_rel_error(n)  = std(objectives{n}.rel_error);
        error_stats.avg_true_error(n) = mean(objectives{n}.true_error);
        error_stats.std_true_error(n) = std(objectives{n}.true_error);
        error_stats.avg_Derror(n)     = mean(objectives{n}.D_error); % Duplicate of avg_error
        error_stats.std_Derror(n)     = std(objectives{n}.D_error);  % Duplicate of std_error
    end
end
