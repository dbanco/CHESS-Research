%% Setup
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\Simulations';
% top_dir = '/cluster/shared/dbanco02/Simulations';

% Simulation name
sim = 'anomaly';
sim_name = [sim,'PoissonNoise4'];
ps_name = 'lineSearchNoise';
mc_dir = 'mcTrials';

load(fullfile(top_dir,sim_name,'IndepISM',[ps_name,'_1']));

% Define poisson dataset
[N,T] = size(B);
[P,K,M,MM] = definePoissonP(N,T);
P.trials = 50;
levels = 0.05:0.05:0.3;
NN = numel(levels);
% [alpha_vals,theta_stds] = genPoissonScales(N,T,levels,sim);
[~,~,theta_stds,~] = genSimDataPoisson(N,T,1,sim);
snrs = (1./levels); snrs(1) = 19;

%% Results of MC Trials
% [awmv_coup,awmv_rse_coup,...
%  awmv_indep,awmv_rse_indep] = awmvMCResults(NN,P.trials,T,theta_stds,sim,...
%                     top_dir,sim_name,mc_dir);

%% Figure4/5b
% AWMV RMSE plots
% fig45b = awmvRmseFigure(snrs,awmv_rse_coup,awmv_rse_indep);


%% Figure 4/5a
% Waterfall plots
nnInd = [2,4,6];
fig45a = waterfallFigure(nnInd,MM,top_dir,sim_name,ps_name);

%% Figure 4/5c
% Show awmvs
[awmv_indep,awmv_coupled] = getAwmvParamSelect(NN,MM,T,top_dir,sim_name,ps_name);

fig45c = awmvParamSelectFigure(awmv_indep,awmv_coupled,theta_stds,snrs);

saveFigures = 0;
if saveFigures
    saveas(fig45a,'fig45a.png','png')
    saveas(fig45b,'fig45b.png','png')
    saveas(fig45c,'fig45c.png','png')
end
