noise_std = 0:0.03:0.30;
n_eta_levels = sqrt(180.*noise_std.^2)/10;
% n_eta_levels = linspace(0.02,0.35,numel(noise_std));

for n_level = 3

    % Parameter selection
    disp('Setup params')
    P.set = 1;

    % Parent directory
    top_dir = 'D:\CHESS_data';
%     top_dir = '/cluster/shared/dbanco02';
    
    % Input dirs
    dset_name = 'gnoise4_nonorm';
    dset_subdir = ['simulated_two_spot_1D_',dset_name,'_',num2str(n_level)];
    indep_name = 'ADMM_Sherman_indep5';
    indep_subdir = [dset_subdir,'_indep1'];
    init_subdir =  [dset_subdir,'_simul_init'];
    
    % Output dirs
    output_name = 'gnoise4_nonorm_coupled_ISM_TVx5';
    output_subdir = [dset_subdir,'_coupled'];
    num_ims = 20;
    
    % Setup directories
    dataset =  fullfile(top_dir,dset_subdir);
    indep_dir = fullfile(top_dir,indep_name,indep_subdir);
    init_dir =  fullfile(top_dir,indep_name,init_subdir);
    output_dir  = fullfile(output_name,output_subdir);

    mkdir([top_dir,output_name])  
    
    % Universal Parameters
    % Ring sampling parameters
    prefix = 'polar_vector';
    baseFileName = 'fista_fit_%i_%i.mat';

    % Load most parameters by loading single output
    load(fullfile(indep_dir,sprintf(baseFileName,1,1)))
    N = numel(P.lambda_values);
    % coupled params
    Pc.initialization = 'simultaneous';
    Pc.preInitialized = 2;
    Pc.rho2 = 10;
    Pc.lambda2 = 0.001;
    Pc.maxIterReg = 1600;
    Pc.tolerance = 1e-10;
    Pc.num_outer_iters = 2;
    Pc.baseFileName = 'fista_fit_%i_%i.mat';
    Pc.num_ims = num_ims;
    Pc.prefix = 'polar_vector';
    Pc.dataset = dataset;

    % Lambda2 values
    lambda2_vals = logspace(-4,-2,30);
    M = numel(lambda2_vals);
    Pc.lambda2_values = lambda2_vals;
    
    % Construct dictionary
    switch P.basis
        case 'norm2'
            A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
    end
    A0 = unshifted_basis_vector_stack_norm2_zpad(P);

    %% Select lambda values
    disp('Selecting lambda values')

    err_select = zeros(N,num_ims);
    l0_select = zeros(N,num_ims);
    l1_select = zeros(N,num_ims);
    for i = 1:N
        fprintf('%i of %i \n',i,N)
        for j = 1:num_ims
            e_data = load(fullfile(indep_dir,sprintf(baseFileName,i,j)),'err','x_hat');
            err_select(i,j) = e_data.err(end);
            l0_select(i,j) = sum(e_data.x_hat(:) > 0);
            l1_select(i,j) = sum(e_data.x_hat(:));
        end
    end
    err_select(err_select > 10^10) = 0;
    l0_select(l0_select > 10^10) = 0;
    l1_select(l1_select > 10^10) = 0;

    % Criterion 
    noise_eta = n_eta_levels(n_level);
    discrep_crit = abs(err_select'-noise_eta);
    [lambda_indices,~] = find(discrep_crit' == min(discrep_crit'));
    param_select = P.lambda_values(lambda_indices);
    Pc.lambda_values = param_select;

    %% Move independent fits to init directory
    mkdir(init_dir)
    for i = 1:num_ims
        src = fullfile(indep_dir,sprintf(baseFileName,lambda_indices(i),i));
        dest = fullfile(init_dir,sprintf(baseFileName,1,i));
        copyfile(src,dest)
    end

    %% Run coupled grid search
    disp('Begin grid search')

    for i = 1:M
        Pc.init_dir = init_dir;
        Pc.output_dirA = [fullfile(top_dir,output_dir),'_',num2str(i),'a'];
        Pc.output_dirB = [fullfile(top_dir,output_dir),'_',num2str(i),'b'];
        Pc.output_dirFinal = [fullfile(top_dir,output_dir),'_',num2str(i),'_','final'];
        mkdir(Pc.init_dir)
        mkdir(Pc.output_dirA)
        mkdir(Pc.output_dirB)
        mkdir(Pc.output_dirFinal)
        Pc.lambda2 = lambda2_vals(i);
        runCoupledISM_TVx_1D(P,Pc)
    end
    
end
