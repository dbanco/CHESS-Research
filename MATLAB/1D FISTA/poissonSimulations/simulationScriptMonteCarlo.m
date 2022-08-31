%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\Simulations';
% top_dir = '/cluster/shared/dbanco02/Simulations';
mkdir(top_dir)

% Simulation name
sim = 'linear';
sim_name = [sim,'PoissonNoise3'];
paramSearch_name = 'lineSearchNoise';
mc_dir = 'mcTrials';

load(fullfile(top_dir,sim_name,'IndepISM',[paramSearch_name,'_1']));

% Define poisson dataset
[N,T] = size(B);
[P,K,M,MM] = definePoissonP(N,T);
P.trials = 100;
levels = 0.05:0.05:0.3;
[alpha_vals,theta_stds] = genPoissonScales(N,T,levels,sim);

%% Run algorithm

% Independent
alg_name = 'IndepISM';
paramSearch_dir  = fullfile(top_dir,sim_name,alg_name);
SimMCPoisson(P,N,K,T,levels,alpha_vals,...
       paramSearch_dir,paramSearch_name,mc_dir,sim)

% Coupled
alg_name = 'CoupledCGTV';
paramSearch_dir  = fullfile(top_dir,sim_name,alg_name); 
SimMCPoissonCoupled(P,N,K,T,levels,alpha_vals,...
    paramSearch_dir,paramSearch_name,mc_dir,sim)

