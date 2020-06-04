noise_std = 0:0.03:0.30;
n_eta_levels = sqrt(180.*noise_std.^2)/100;
% n_eta_levels = linspace(0.02,0.35,numel(noise_std));

for n_level = 3

    % Parameter selection
    disp('Setup parms')
    P.set = 1;

    dset_name = 'gnoise4_nonorm';
    num_ims = 20;
    
%     datadir = '/cluster/shared/dbanco02/';
%     dataset = ['/cluster/home/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'/'];
%     base_dir = '/cluster/shared/dbanco02/ADMM_Sherman_indep3/';
%     indep_dir = [base_dir,'simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_indep/'];
%     init_dir =  [base_dir,'simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_simul_init'];
%     output_dir = ['gnoise4_coupled_ISM1/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_coupled'];
    
    datadir = 'D:\CHESS_data\';
    dataset =  [datadir,'simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'\'];
    base_dir = [datadir,'ADMM_Sherman_indep3\'];
    indep_dir = [base_dir,'simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_indep1\'];
    init_dir =  [base_dir,'simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_simul_init'];
    output_dir  = ['gnoise4_coupled_ISM1\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_coupled'];

    mkdir([datadir,'gnoise4_coupled_ISM1'])    
    % Universal Parameters
    % Ring sampling parameters
    prefix = 'polar_vector';
    baseFileName = 'fista_fit_%i_%i.mat';

    % Load most parameters by loading single output
    load([indep_dir,sprintf(baseFileName,1,1)])
    N = numel(P.lambda_values);
    % coupled params
    Pc.initialization = 'simultaneous';
    Pc.preInitialized = 2;
    Pc.rho2 = 1;
    Pc.lambda2 = 0.001;
    Pc.maxIterReg = 800;
    Pc.num_outer_iters = 10;
    Pc.baseFileName = 'fista_fit_%i_%i.mat';
    Pc.num_ims = num_ims;
    Pc.prefix = 'polar_vector';
    Pc.dataset = dataset;
    Pc.distScale = 0;

    % Lambda2 values
    lambda2_vals = logspace(-2,-1,10);
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
            err_select(i,j) = e_data.err(end-1)*norm(polar_image);
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
        Pc.output_dirA = [datadir,output_dir,'_',num2str(i),'a'];
        Pc.output_dirB = [datadir,output_dir,'_',num2str(i),'b'];
        mkdir(Pc.init_dir)
        mkdir(Pc.output_dirA)
        mkdir(Pc.output_dirB)
        Pc.gamma = gamma_vals(i);
        runCoupledISM_TVx_1D(P,Pc)
    end
    
end
