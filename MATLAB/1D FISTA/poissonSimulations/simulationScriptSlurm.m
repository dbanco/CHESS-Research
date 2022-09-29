%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'D:\Simulations';
top_dir = '/cluster/shared/dbanco02/Simulations';
mkdir(top_dir)

% Simulation name
sim = 'anomaly';
sim_name = [sim,'PoissonNoise6'];
file_name = 'lineSearchNoise';

% Define poisson dataset
N = 101;
T = 30;
[P,K,M,MM] = definePoissonP(N,T,sim);
snr_levels = [10,5,2.5,2,1.75,1.5,1.25,1.1];
[alpha_vals,theta_stds] = genPoissonScales2(N,T,snr_levels,sim);
P.alpha_vals = alpha_vals;
% [Bn,B,theta_stds,rel_err] = genSimDataPoisson2(N,T,alpha_vals(4),sim);

%% Run algorithm

% Independent
alg_name = 'IndepISM';
indep_dir  = fullfile(top_dir,sim_name,alg_name);
SimParamSearchPoisson2(P,M,snr_levels,alpha_vals,...
                       indep_dir,file_name,sim)


% Coupled
alg_name = 'CoupledCGTV';
coupled_dir  = fullfile(top_dir,sim_name,alg_name);       
SimParamSearchCoupledPoissonSlurm(P,MM,snr_levels,...
                       coupled_dir,indep_dir,file_name,sim_name)


