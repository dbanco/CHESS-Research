% Compute metrics 
sigmas = 0:0.01:0.1;
dataset = "dissertation_adjust2";
% Local OF recovery
load('E:\MCDLOF_processing\Outputs_8_19_local_dissertation_adjust2_log_Dmcdl0_Xzeros0_recenter0\results_trial_1_sig_2\output_j41_10_5_sig_1.00e-02_lam1_2.21e-02_lam2_2.08e-01_lam3_2.00e-03.mat')

opt = outputs.opt;
lambda2 = outputs.lambda2;
lambda3 = opt.Smoothness;
HSiters = opt.HSiters;

% --- Compute metrics of OF run ---
[error_of, rel_error_of, of_penalty_of, hs_penalty_of,D_error_of] = ...
    compute_recov_metrics(outputs, dataset, sigmas(2), lambda2, lambda3, HSiters);

% --- Compute metrics for independent run ---
load("E:\MCDLOF_processing\Outputs_8_18_trials_indep_dissertation_adjust2_log_Dflat0_Xzeros0_recenter0\results_trial_1_sig_2\output_j41_1_1_sig_1.00e-02_lam1_2.21e-02_lam2_0.00e+00_lam3_0.00e+00.mat")
[error_indep, rel_error_indep, of_penalty_indep, hs_penalty_indep,D_error_indep] = ...
    compute_recov_metrics(outputs, dataset, sigmas(2), lambda2, lambda3, HSiters);

% --- Organize results ---
labels = {'Error', 'Relative Error', 'OF Penalty', 'HS Penalty', ...
          'D Error'};

indep_vals = [error_indep, rel_error_indep, of_penalty_indep, hs_penalty_indep,D_error_indep];
of_vals = [error_of, rel_error_of, of_penalty_of, hs_penalty_of,D_error_of];

% --- Display neatly ---
fprintf('\n%-15s | %-12s | %-12s\n', 'Metric', 'Independent', 'OF');
fprintf('%s\n', repmat('-', 1, 45));

for i = 1:numel(labels)
    fprintf('%-15s | %-12.4f | %-12.4f\n', labels{i}, indep_vals(i), of_vals(i));
end

[~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigmas(2),dataset);