%% Setup
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\Simulations';
% top_dir = '/cluster/shared/dbanco02/Simulations';

% Simulation name
sim = 'anomaly';
sim_name = [sim,'PoissonNoise3'];
ps_name = 'lineSearchNoise';
mc_dir = 'mcTrials';

load(fullfile(top_dir,sim_name,'IndepISM',[ps_name,'_1']));

% Define poisson dataset
[N,T] = size(B);
[P,K,M,MM] = definePoissonP(N,T);
P.trials = 50;
levels = 0.05:0.05:0.3;
% [alpha_vals,theta_stds] = genPoissonScales(N,T,levels,sim);
[~,~,theta_stds,~] = genSimDataPoisson(N,T,1,sim);


%% Results of MC Trials
[awmv_coup,awmv_rse_coup,...
 awmv_indep,awmv_rse_indep] = awmvMCResults(NN,nTrials,T,theta_stds,sim);

%% Figure4/5b
% AWMV RMSE plots
fig45b = awmvRmseFigure(levels,awmv_rse_coup,awmv_rse_indep);


%% Figure 4/5a
% Waterfall plots
nnInd = [2,4,6];
fig45a = waterfallFigure('linear',nnInd);

%% Figure 4/5c
% Show awmvs
[awmv_indep,awmv_coupled] = getAwmvParamSelect(sim,NN,T);
fig45c = awmvParamSelectFigure(awmv_indep,awmv_coupled,theta_stds,levels);

