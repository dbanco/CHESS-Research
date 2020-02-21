%% Parameter selection
clear all
close all
disp('Setup parms')
P.set = 1;
datadir = '/cluster/shared/dbanco02/';
dataset = '/cluster/home/dbanco02/mmpad_polar/ring1_zero/';
indep_dir = '/cluster/shared/dbanco02/mmpad_1D_indep_param_8/';
output_dir = '/cluster/shared/dbanco02/mmpad_1D_coupled_param_8/';
mkdir(output_dir)

num_ims = 500;
prefix = 'mmpad_img';
baseFileName = 'fista_fit_%i_%i.mat';

% Lambda values
lambda_vals = logspace(-3,1,30); 
N = numel(lambda_vals);

% Load most parameters by loading single output
load([indep_dir,sprintf(baseFileName,1,1)])

% coupled params
Pc.initialization = 'simultaneous';
Pc.preInitialized = 2;
Pc.wLam = 10;
Pc.gamma = 1;
Pc.maxIterReg = 800;
Pc.num_outer_iters = 10;
Pc.baseFileName = 'fista_fit_%i_%i.mat';
Pc.num_ims = num_ims;
Pc.prefix = 'mmpad_img';
Pc.dataset = [dataset];

% Gamma values
gamma_vals = [ 0.01,0.025 0.05,0.075,0.1,0.15,0.2,0.0075,0.005,0.0025,0.001,0.00075,0.0005]; 
M = numel(gamma_vals);

%% Select lambda values
disp('Selecting lambda values')

err_select = zeros(N,num_ims);
l0_select = zeros(N,num_ims);
l1_select = zeros(N,num_ims);
for i = 1:N
    fprintf('%i of %i \n',i,N)
    for j = 1:num_ims
        e_data = load(fullfile(indep_dir,sprintf(baseFileName,i,j)),'err','x_hat');
        err_select(i,j) = e_data.err(end-1);
        l0_select(i,j) = sum(e_data.x_hat(:) > 0);
        l1_select(i,j) = sum(e_data.x_hat(:));
    end
end
err_select(err_select > 10^10) = 0;
l0_select(l0_select > 10^10) = 0;
l1_select(l1_select > 10^10) = 0;

% Criterion 
noise_eta = 0.10;
discrep_crit = abs(err_select'-noise_eta);

[lambda_indices,~] = find(discrep_crit' == min(discrep_crit'));
param_select = lambda_vals(lambda_indices);
Pc.lambda_values = param_select;

%% Move independent fits to init directory (done)

init_dir = [datadir,'mmpad_1D_coupled_simul_init'];
% mkdir(init_dir)
% for i = 1:num_ims
%     src = fullfile(indep_dir,sprintf(baseFileName,lambda_indices(i),i));
%     dest = fullfile(init_dir,sprintf(baseFileName,1,i));
%     copyfile(src,dest)
% end

%% Run coupled grid search

disp('Begin grid search')

for i = 11:13
    Pc.init_dir = init_dir;
    Pc.output_dirA = [datadir,'mmpad_1D_coupled_simul_',num2str(i),'a'];
    Pc.output_dirB = [datadir,'mmpad_1D_coupled_simul_',num2str(i),'b'];
    mkdir(Pc.init_dir)
    mkdir(Pc.output_dirA)
    mkdir(Pc.output_dirB)
    Pc.gamma = gamma_vals(i);
    
    runCoupledFISTA_1D(P,Pc)
end
