%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = logspace(-2,0,150);
lambdaOFVals = [0 logspace(-4,0,20)];
lambdaHSVals = [0 logspace(-4,0,10)];

selected_lam_s = [41,41,66,84,92,96];

% selected_lam_of = [1,1,1,1,1,1];
% selected_lam_hs = [1,1,1,1,1,1];

% log-penalty based parameter selection
% selected_lam_of = [21,21,10,9,11,12];
% selected_lam_hs = [3,3,6,7,10,6];

%of+log penalty based selection
% selected_lam_of = [1,21,21,15,21,19]; 
% selected_lam_hs = [1,3,2,3,4,7];

%fixed parameters
selected_lam_of = [1,10,10,10,10,10,10]; 
selected_lam_hs = [1,6,6,6,6,6,6];

% Experiment Setup
sigmas = 0:0.01:0.1;
NN = numel(sigmas);

fig_num = 22;

datasets = {'dissertation_adjust2'};
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

test_name = ['Outputs_8_27_test_of_trials_',dataset,'_',opt.Penalty,...
    '_Dmcdl',num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

test_name2 = ['Outputs_9_4_of_trials_',dataset,'_',opt.Penalty,...
    '_Dmcdl',num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

topDir = ['E:\MCDLOF_processing\',test_name];
topDir2 = ['E:\MCDLOF_processing\',test_name2];

useMin = 1;
sig_ind = 2:6;
num_trials = 50;

[objectives_of,objectives_indep] = compute_objectives_trials(sig_ind,sigmas,...
    dataset,useMin,num_trials,topDir,topDir2,...
    selected_lam_s,selected_lam_of,selected_lam_hs,...
    lambdaVals,lambdaOFVals,lambdaHSVals);

%% Plot Results
fig_dir = 'C:\Users\dpqb1\Documents\MCDL Paper';
data_name = 'sim2';

[meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);

error_stats_indep = compute_error_stats(objectives_indep,sig_ind);
error_stats_of = compute_error_stats(objectives_of,sig_ind);

trials_figures_subplot(fig_dir,num_trials,data_name,meanSNR,error_stats_indep,error_stats_of)