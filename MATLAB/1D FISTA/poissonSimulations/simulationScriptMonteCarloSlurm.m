%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\Simulations';
% top_dir = '/cluster/shared/dbanco02/Simulations';
mkdir(top_dir)

% Simulation name
sim = 'anomaly';
sim_name = [sim,'PoissonNoise6'];
ps_name = 'lineSearchNoise';
mc_dir = 'mcTrials';

load(fullfile(top_dir,sim_name,'IndepISM',[ps_name,'_1']));

% Define poisson dataset
[N,T] = size(B);
snr_levels = [10,5,2.5,2,1.75,1.5,1.25,1.1];
P.trials = 100;
alpha_vals = P.alpha_vals;
theta_stds = P.theta_stds;

%% Run algorithm

% Independent
alg_name = 'IndepISM';
indepPS_dir  = fullfile(top_dir,sim_name,alg_name);
SimMCPoisson2(P,snr_levels,alpha_vals,...
              indepPS_dir,ps_name,mc_dir,sim)

% Coupled
alg_name = 'CoupledCGTV';
coupledPS_dir = fullfile(top_dir,sim_name,alg_name); 
SimMCPoissonCoupledSlurm(snr_levels,indepPS_dir,coupledPS_dir,ps_name,mc_dir,sim_name)

