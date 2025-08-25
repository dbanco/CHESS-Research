% Compute metrics 
sigmas = 0:0.01:0.1;
dataset = "dissertation_adjust2";
% Local OF recovery
load("E:\MCDLOF_processing\Outputs_8_25_of_dissertation_adjust2_log_Dmcdl0_Xzeros0_recenter0\results_trial_1_sig_2\output_j39_21_3_sig_1.00e-02_lam1_3.24e-02_lam2_1.00e+00_lam3_2.78e-04.mat")

opt = outputs.opt;
lambda1 = outputs.lambda;
lambda2 = outputs.lambda2;
lambda3 = opt.Smoothness;
HSiters = opt.HSiters;

% --- Compute metrics of OF run ---
[outVector_of,outLabels] = ...
    compute_recov_metrics(outputs, dataset, sigmas(2),lambda1, lambda2, lambda3, HSiters);

% --- Compute metrics for independent run ---
load("E:\MCDLOF_processing\Outputs_8_25_indep_dissertation_adjust2_log_Dflat0_Xzeros0_recenter0\results_trial_1_sig_2\output_j39_1_1_sig_1.00e-02_lam1_3.24e-02_lam2_0.00e+00_lam3_0.00e+00.mat")
[outVector_indep,outLabels] = ...
    compute_recov_metrics(outputs, dataset, sigmas(2), lambda1, lambda2, lambda3, HSiters);

% --- Display neatly ---
fprintf('\n%-15s | %-12s | %-12s\n', 'Metric', 'Independent', 'OF');
fprintf('%s\n', repmat('-', 1, 45));

for i = 1:numel(outLabels)
    fprintf('%-15s | %-12.4f | %-12.4f\n', outLabels{i}, outVector_indep(i), outVector_of(i));
end

[~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigmas(2),dataset);