%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\Simulations';
% top_dir = '/cluster/shared/dbanco02';
mkdir(top_dir)

% Simulation name
sim_name = 'linearPoissonNoise';
file_name = 'lineSearchNoise';

% Define poisson dataset
N = 101;
T = 30;
[P,K,M,MM] = definePoissonP(N,T);
levels = 0.05:0.05:0.3;
alpha_vals = genPoissonScales(N,T,levels,'linear');

% Independent
alg_name = 'IndepISM';
indep_dir  = fullfile(top_dir,sim_name,alg_name);
linearSimParamSearchPoisson(P,N,M,K,T,levels,alpha_vals,...
                                indep_dir,file_name)

% Coupled
alg_name = 'CoupledCGTV';
coupled_dir  = fullfile(top_dir,sim_name,alg_name);       
linearSimParamSearchCoupledPoisson(P,N,MM,K,T,levels,alpha_vals,...
                                coupled_dir,indep_dir,file_name)
                            