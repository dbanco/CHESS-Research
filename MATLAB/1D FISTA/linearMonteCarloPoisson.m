%% Setup
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'D:\CHESS_data\';
top_dir = '/cluster/home/dbanco02/data';

num_ims = 50;
N = 101;
K = 20;
M = 50;
NN = 1;
T = num_ims;
zPad = 0;
zMask = [];


theta_stds1 = linspace(1,15,T);

%% Parameter Selection Independent
close all
lambda_select = zeros(NN,T);

for nn = 1:NN
    dset_name = ['singlePeak_noise_Poisson'];
    indep_name = '_indep_ISM';
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    indep_subdir = [dset_name,indep_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    indep_dir  = fullfile(top_dir,indep_subdir);
    load(fullfile(indep_dir,[dset_name,'_',num2str(num_ims),'_','all']))
    if( nn==1 )
        % Construct dictionary
        A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
    end
    
    % Compute outputs
    rse_indep = zeros(M,T);
    mse_indep = zeros(M,T);
    l1_norm = zeros(M,T);
    for i = 1:M
        for time = 1:T
            x = X_indep(:,:,i,time);
            fit = Ax_ft_1D(A0ft_stack,x);
            l1_norm(i,time) = sum(abs(x(:)));
            mse_indep(i,time) = norm(fit-B(:,time));
            rse_indep(i,time) = norm(fit-B(:,time))/norm(B(:,time));   
        end
    end

    % Parameter Selection
    select_ind = zeros(T,1);
    for time = 1:T
        crit = abs(l1_norm(:,time)*0.5).^2 + abs(mse_indep(:,time)).^2;
        select_ind(time) = find( (crit == min(crit)),1 );
        lambda_select(nn,time) = P.lambda_values(select_ind(time));
    end
end
   

for nn = 1:NN
    % Input dirs
    dset_name = ['singlePeak_noise_MC_poisson'];

    % Output dirs
    output_name = '_indep_ISM';
    output_subdir = [dset_name,output_name];

    % Setup directories
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    mkdir(output_dir)

    num_ims = 50;

    % File Parameters
    P.baseFileName = 'indep_fit_%i_%i.mat';
    P.dataset = dataset;

    % Data/Dictionary Parameters
    % Zero padding and mask

    N = 101;
    K = 20;
    M = 50;
    T = num_ims;
    zPad = 0;
    zMask = [];

    P.dataScale = 1;
    P.lambda_values = logspace(-3,1,M);
    P.num_theta = N;
    P.sampleDims = [T,1];
    P.num_ims = T;
    P.basis = 'norm2';
    P.cost = 'l1';
    P.num_var_t = K;
    P.var_theta = linspace(0.5,25,P.num_var_t).^2;

    % algorithm parameters

    P.params.rho1 = 1;
    % P.params.lambda1 = 0.0001;
    P.params.tau = 1.05;
    P.params.mu = 2;
    P.params.adaptRho = 1;
    P.params.alpha = 1.8;
    P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
    P.params.maxIter = 50;
    P.params.tolerance = 1e-8;
    P.params.isNonnegative = 1;
    P.params.zeroPad = zPad;
    P.params.zeroMask = zMask;
    P.params.plotProgress = 0;
    P.params.verbose = 1;

    P.params.conjGradIter = 100;
    P.params.tolerance = 1e-8;
    P.params.cgEpsilon = 1e-3;

    % Construct dictionary
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

    % Define data
    theta_stds1 = linspace(1,15,T);
    

    %% Independent Solution
    trials = 100;
    x_init = zeros(N,K);
    X_indep = zeros(N,K,trials,T);
    B = zeros(N,T,trials);
    % Run 100 trials
    for i = 1:trials
        P.set = i;
        % Regenerate noisy data
        for t = 1:T
            b = gaussian_basis_1D( N, N/2, theta_stds1(t)^2);
            rms = sqrt(sum(b.^2)/N);
            b = b*100; % + randn(N,1)*noise_std(nn);
            bn = poissrnd(b);
            B(:,t,i) = bn;
            % Solve
            P.params.lambda1 = lambda_select(nn,t);
            [x_hat,obj,err,l1_norm,~] = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_init,P.params);  
            X_indep(:,:,i,t) = x_hat;
        end
    end
    save(fullfile(output_dir,[dset_name,'_',num2str(nn),'_','all']),...
        'B','X_indep','P');
end
