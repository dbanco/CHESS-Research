%% Parameter selection
clear all
close all
P.set = 1;
datadir = 'D:\CHESS_data\';
for d_num = 1:8
    dataset_num = num2str(d_num);
    dataset = ['D:\CHESS_data\simulated_two_spot_1D_noise2_',dataset_num,'\'];
    param_dir = 'D:\CHESS_data\param_search_1D\';
    num_ims = 10;

    % Universal Parameters
    % Ring sampling parameters
    prefix = 'polar_vector';
    load([dataset,prefix,'_1.mat']);
    P.num_theta = size(polar_vector,1);
    P.dtheta = 1;
    P.sampleDims = [num_ims,1];

    % Basis function variance parameters
    P.basis = 'norm2';
    P.cost = 'l1';
    P.num_var_t = 15;
    P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;

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
    load([param_dir,'lambda_select_noise8'])
    Pc.lambda_values = param_select;
    Pc.initialization = 'simultaneous';
    Pc.wLam = 25;
    Pc.gamma = 1;
    Pc.maxIterReg = 800;
    Pc.num_outer_iters = 10;
    Pc.baseFileName = 'fista_fit_%i_%i.mat';
    Pc.num_ims = num_ims;
    Pc.prefix = 'polar_vector';
    Pc.dataset = [dataset];
    % Gamma values
    gamma_vals = [ 0.05,0.1, 0.2,0.01, 0.025, 0.035]; 
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
    %% Run grid search
    for i = 6
        Pc.init_dir = [datadir,'noise_1D_coupled_simul_',dataset_num,'_init',num2str(i)];
        Pc.output_dirA = [datadir,'noise_1D_coupled_simul_',dataset_num,'_',num2str(i),'a'];
        Pc.output_dirB = [datadir,'noise_1D_coupled_simul_',dataset_num,'_',num2str(i),'b'];
        mkdir(Pc.init_dir)
        mkdir(Pc.output_dirA)
        mkdir(Pc.output_dirB)
        Pc.gamma = gamma_vals(i);
        runCoupledFISTA_1D(P,Pc)
    end
end
% save([param_dir,'param_search_coupled_discrep_8.mat'],'err_gcv','obj_gcv','P','lambda_vals')
