%% Setup
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'D:\Simulations';
top_dir = '/cluster/shared/dbanco02/Simulations';

% Simulation name
sim = 'linear';
sim_name = [sim,'PoissonNoise6'];
ps_name = 'lineSearchNoise';
mc_dir = 'mcTrials';

load(fullfile(top_dir,sim_name,'IndepISM',[ps_name,'_1']));
load(fullfile(top_dir,sim_name,'IndepISM','mcTrials_1','trial_1'));
% Define poisson dataset
[N,T] = size(B);
[P,K,M,MM] = definePoissonP(N,T,sim);
P.trials = 2;
snr_levels = [10,5,2.5,2,1.75,1.5,1.25,1.1];
NN = numel(snr_levels);

theta_stds = P.theta_stds;
snrs = snr_levels;

                 

%% Figure 4/5a
% Waterfall plots
nnInd = [2,5,7];
fig45a = waterfallFigure(nnInd,MM,top_dir,sim_name,ps_name);

%% Figure 4/5c
% Show awmvs
[awmv_indep,awmv_coupled] = getAwmvParamSelect(NN,MM,T,top_dir,sim_name,ps_name);

fig45c = awmvParamSelectFigure(awmv_indep,awmv_coupled,theta_stds,snr_levels);

%% Results of MC Trials
[awmv_coup,awmv_rse_coup,...
 awmv_indep,awmv_rse_indep] = awmvMCResults(NN,2,T,theta_stds,sim,...
                                              top_dir,sim_name,mc_dir);

%% Figure4/5b
% AWMV RMSE plots
fig45b = awmvRmseFigure(snrs,awmv_rse_coup,awmv_rse_indep);

saveFigures = 0;
if saveFigures
    saveas(fig45a,'fig45a.png','png')
    saveas(fig45b,'fig45b.png','png')
    saveas(fig45c,'fig45c.png','png')
end
