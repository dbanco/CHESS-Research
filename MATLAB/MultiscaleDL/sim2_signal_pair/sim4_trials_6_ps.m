%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = logspace(-2,0,150);
lambdaOFVals = [0 logspace(-4,0,20)];
lambdaHSVals = [0 logspace(-4,0,10)];

selected_lam_s = [66,88,115,123,132];
selected_lam_of = [7,7,4,8,7];
selected_lam_hs = [4,11,11,7,11];

lambda_eval = [9,6];

% Experiment Setup
% Noise level
dataset = 'gaus_loren_nooverlap';
sig_ind = 1:5;
SNRs = [20,16,12,8,4];
sigmas = zeros(numel(SNRs),1);
for i = sig_ind
    sigmas(i) = SNRtoSigma(SNRs(i),dataset);
end

NN = numel(sigmas);
fig_num = 22;

datasets = {'gaus_loren_nooverlap'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true','mcdl'};
dfixes = {0,1};
recenters = {0,1};

% Setup Dataset
s0 = 1;
s1 = 2;
s2 = 1;
s3 = 2;
s4 = 1;
s5 = 1;
dataset = datasets{s0};
opt.Penalty = penalties{s1};
opt.coefInit = xinits{s2};
opt.dictInit = dinits{s3};
opt.Dfixed = dfixes{s4};
opt.Recenter = recenters{s5};
opt.Xfixed = 0;

test_name = ['Outputs_9_8_indep_',dataset,'_',opt.Penalty,...
    '_Dflat',num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

test_name2 = ['Outputs_9_8_of_',dataset,'_',opt.Penalty,...
    '_Dflat',num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

topDir = ['E:\MCDLOF_processing\',test_name];
topDir2 = ['E:\MCDLOF_processing\',test_name2];

useMin = 1;
num_trials = 1;

[objectives_of,objectives_indep] = compute_objectives_trials(sig_ind,sigmas,...
    dataset,useMin,num_trials,topDir,topDir2,...
    selected_lam_s,selected_lam_of,selected_lam_hs,...
    lambdaVals,lambdaOFVals,lambdaHSVals,lambda_eval);

%% Plot Results
fig_dir = 'C:\Users\dpqb1\Documents\MCDL Paper';
data_name = 'sim4';

error_stats_indep = compute_error_stats(objectives_indep,sig_ind);
error_stats_of = compute_error_stats(objectives_of,sig_ind);

trials_figures_subplot_sim2(fig_dir,num_trials,data_name,SNRs,error_stats_indep,error_stats_of)