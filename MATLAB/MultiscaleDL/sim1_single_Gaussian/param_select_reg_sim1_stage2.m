lambdaVals = logspace(-1.5,1,150);
lambdaRegVals = [0 logspace(-2,2,99)];

fig_num = 22;

datasets = {'pseudo-voigt_unmatched'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true','mcdl'};
dfixes = {0,1};
recenters = {0,1};

% Noise level
dataset = datasets{1};
sig_ind = 1:5;
SNRs = [20,17,14,11,8];
sigmas = zeros(numel(SNRs),1);
for i = sig_ind
    sigmas(i) = SNRtoSigma(SNRs(i),dataset);
end
NN = numel(sigmas);

% Setup Dataset
s = [1,2,1,2,1,1];
opt.Penalty = penalties{s(2)};
opt.coefInit = xinits{s(3)};
opt.dictInit = dinits{s(4)};
opt.Dfixed = dfixes{s(5)};
opt.Recenter = recenters{s(6)};
opt.Xfixed = 0;

prefix = ['10_14_reg_softmin_LBFGS_',dataset];
topDir = ['E:\MCDLOF_processing\Outputs_',prefix,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed)];

% criterion = 'discrepancy';
% criterion = 'truth_error';
% criterion = 'relaxed discrepancy';
% criterion = 'discrepancy range';
% criterion = 'score';
% criterion = 'l-curve-MDF3b';
% criterion = 'wass_dist';
criterion = 'x_metric2';
% criterion = 'discrepancy range of-hs-log';


selected_lam_s_vec = zeros(NN,1);
selected_lam_reg_vec = zeros(NN,1);
selected_inds = zeros(NN,1);
objectives = cell(NN,1);

useMin = 1;
relax_param = 1.1; % for discrepancy range
makeFigures = true;

for n = sig_ind
    inDir = [topDir,'\results_trial_1_sig_',num2str(n)];
    resultFile = [topDir,'\results_sig_',num2str(n),'.mat'];
    % if exist(resultFile,'file')
    %     load(resultFile)
    % else
    results = compute_metrics_reg(inDir,sigmas(n),dataset,useMin);
        % save(resultFile,'results')
    % end
    [lambda_all,objective] = param_select_mcdl_regularized(results,criterion,sigmas(n),dataset,relax_param,fig_num);
    selected_lam_all_vec(n,:) = lambda_all;
    selected_inds(n,1) = find(selected_lam_all_vec(n,1) == lambdaVals);
    selected_inds(n,2) = find(selected_lam_all_vec(n,2) == lambdaRegVals);

    j_s = selected_inds(n,1);
    j_reg = selected_inds(n,2);

    if makeFigures
        outputs = loadOutputFile(inDir,selected_inds(n,:));
        suffix = sprintf('_j%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e',...
                      j_s,j_reg,sigmas(n),outputs.lambda,outputs.lambda2);
        psFigDir = fullfile(topDir,['ps_',criterion,'_',num2str(relax_param),...
                                          '_useMin',num2str(useMin)]);
        generateFiguresToy1zpad_center(psFigDir,outputs,suffix,[4,8],useMin);
    end
    objectives{n} = objective;
end

% Plot Evaluation Metrics
prefix2 = ['10_10_softmin_LBFGS_',dataset];
topDirIndep = ['E:\MCDLOF_processing\Outputs_',prefix2,'_',opt.Penalty,'_Dflat0_Xzeros0'];
topDirReg = ['E:\MCDLOF_processing\Outputs_',prefix,'_',opt.Penalty,'_Dflat0_Xzeros0'];
num_trials = 1;
[objectives_of,objectives_indep] = compute_objectives_metrics(sig_ind,sigmas,...
    dataset,topDirIndep,topDirReg,selected_inds(:,1),selected_inds(:,2));

data_name = 'sim2';

error_stats_indep = compute_error_stats2(objectives_indep,sig_ind);
error_stats_of = compute_error_stats2(objectives_of,sig_ind);

plot_eval_metrics(psFigDir,num_trials,data_name,SNRs,error_stats_indep,error_stats_of)
