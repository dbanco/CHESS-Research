%% Parameter selection
% Run grid search
clear all
close all
P.set = 1;
datadir = 'D:\CHESS_data\';


dataset = ['D:\MMPAD_data\ring1_zero_subset\'];
param_dir = 'D:\CHESS_data\param_search_mmpad\';
num_ims = 10;

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_image';
load([dataset,prefix,'_1.mat']);
P.num_theta = size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r).^2;
% Zero padding and mask\
zPad = [0,0];
zMask = [];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1;
params.t_k = 1;
params.lambda = 0;
params.wLam = 25;
params.gamma = 1;
params.maxIterReg = 800;
params.beta = 1.2;
params.maxIter = 800;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;

P.params = params;

% coupled params
load([param_dir,'lambda_select_mmpad_subset1'])
Pc.lambda_values = param_select;
Pc.initialization = 'simultaneous';
Pc.wLam = 25;
Pc.gamma = 1;
Pc.maxIterReg = 800;
Pc.num_outer_iters = 10;
Pc.baseFileName = 'fista_fit_%i_%i.mat';
Pc.num_ims = num_ims;
Pc.prefix = 'polar_image';
Pc.dataset = [dataset];

% Gamma values
gamma_vals = [ 0.01,0.02,0.03 0.04 0.05,0.08,0.1,0.15,0.2]; 
N = numel(gamma_vals);

% Select noise level
obj_gcv = zeros(N,num_ims);
err_gcv = zeros(N,num_ims);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
end
A0 = unshifted_basis_vector_stack_norm2(P);

for i = 1:N
    Pc.init_dir = [datadir,'mmpad_subset_coupled_simul_',dataset_num,'_init',num2str(i)];
    Pc.output_dirA = [datadir,'mmpad_subset_coupled_simul_',dataset_num,'_',num2str(i),'a'];
    Pc.output_dirB = [datadir,'mmpad_subset_coupled_simul_',dataset_num,'_',num2str(i),'b'];
    mkdir(Pc.init_dir)
    mkdir(Pc.output_dirA)
    mkdir(Pc.output_dirB)
    Pc.gamma = gamma_vals(i);
    runCoupledFISTA(P,Pc)
end

%% Get fit quantities
clear all
close all
P.set = 1;
datadir = 'D:\CHESS_data\';

% Gamma values
gamma_vals = [ 0.01,0.02,0.03 0.04 0.05,0.08,0.1,0.15,0.2]; 
N = numel(gamma_vals);

obj_values = zeros(8,N);
err_values = zeros(8,N);
err_est_values = zeros(8,1);
norm_values = zeros(8,1);

for i = 1:N
    norm_total = 0;
    err_total = 0;
    est_err_total = 0;
    obj_total = 0;
    for j = 1:10
        dataset_num = num2str(d_num);
        result_dir = [datadir,'noise_1D_coupled_simul_',dataset_num,'_',num2str(i),'a'];
        obj_name = ['objective_10_',num2str(j),'.mat'];
        load(fullfile(result_dir,obj_name))
        obj_total = obj_total + obj(end);

        load(fullfile(result_dir,sprintf('fista_fit_1_%i.mat',j)))
        err_total = err_total + err(end).^2

        % Estimate error
        b = polar_image;
        bn = b/norm(b(:));
        kernel = [0.0001 0.0003 0.0012 0.0042 0.0127 0.0329,...
                  0.0740 0.1434 0.2399 0.3466 0.4322 0.4652];
        kernel = [kernel,fliplr(kernel(1:end-1))];
        kernel = kernel./sum(kernel(:));
        bn_hat = conv(bn,kernel,'same');

        norm_total = norm_total + norm(b(:))^2;
        est_err_total = est_err_total + norm(bn-bn_hat)^2;
    end
    norm_values(d_num) = sqrt(norm_total);
    err_est_values(d_num) = sqrt(est_err_total);
    obj_values(d_num,i) = obj_total;
    err_values(d_num,i) = sqrt(err_total);
end

%% Plot out L-curves
figure(1)
for d_num = 1:8
    plot(obj_values(d_num,:),'o-')
    hold on
end
ylabel('Objective Values')
legend('1','2','3','4','5','6','7','8','9')

figure(2)
for d_num = 1:8
    plot(err_values(d_num,:),'o-')
    hold on
end
ylabel('Error')
xlabel('\lambda_W')
legend('1','2','3','4','5','6','7','8','9')
%% Ty different
