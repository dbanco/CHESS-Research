%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'D:\CHESS_data\';
top_dir = '/cluster/home/dbanco02/data';


num_ims = 30;
N = 101;
K = 20;
M = 50;
T = num_ims;
zPad = 0;
zMask = [];

theta_stds1 = [7*ones(1,T/2),12*ones(1,T/2)]';

%% Paramter Selection Independent
close all
NN = 1;
lambda_select = zeros(NN,T);
dict_init = 1;
for nn = 1:NN
    % Load Independent and compute awmv
    dset_name = ['anomaly_noise_Poisson'];
    indep_name = '_indep_ISM';
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    indep_subdir = [dset_name,indep_name];
    indep_dir  = fullfile(top_dir,indep_subdir);
    load(fullfile(indep_dir,[dset_name,'_',num2str(M),'_','all']))
    if( dict_init )
        % Construct dictionary
        A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
        dict_init = 0;
    end
    mse_indep = zeros(M,T);
    l1_norm = zeros(M,T);
    for i = 1:M
        for time = 1:T
            x = X_indep(:,:,i,time);
            fit = Ax_ft_1D(A0ft_stack,x);
            l1_norm(i,time) = sum(abs(x(:)));
            mse_indep(i,time) = norm(fit-B(:,time));  
        end
    end

    % Parameter Selection
    select_ind = zeros(T,1);
    
    for time = 1:T
        crit = abs(l1_norm(:,time)*0.60).^2 + abs(mse_indep(:,time)).^2;
        select_ind(time) = find( (crit == min(crit)),1 );
        lambda_select(nn,time) = P.lambda_values(select_ind(time));
    end
end

%% Paramter Selection Coupled
MM = 20;
gamma_values = [logspace(-2.5,-1,MM);
              logspace(-2.2,-0.8,MM);
              logspace(-2.2,-0.8,MM);
              logspace(-2.2,-0.5,MM);
              logspace(-1.5,-0.5,MM);
              logspace(-1.5,-0.5,MM);
              logspace(-1.5,-0.5,MM);
              logspace(-1.5,-0.3,MM);
              logspace(-1.5,-0.3,MM);
              logspace(-1.5,-0.3,MM);
              logspace(-1.5,-0.3,MM)];
          
awmv_all = zeros(MM,T,NN);
awmv_rmse = zeros(MM,NN);
mse = zeros(MM,NN);
tv_penalty = zeros(MM,NN);
dict_init = 1;

for nn = 1:NN
    % Load coupled solution
    dset_name = ['anomaly_noise_Poisson'];
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    output_dir  = fullfile(top_dir,output_subdir);
    for i = 1:MM
        load(fullfile(output_dir,[dset_name,'_',num2str(i),'_','time',num2str(1)]))
        if(dict_init)
            % Construct dictionary
            A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
            dict_init = 0;
        end
        for time = 1:T
            x = X_hat(:,:,time);
            fit = Ax_ft_1D(A0ft_stack,x);
            mse(i,nn) = mse(i,nn) + norm( B(:,time)-fit ).^2;
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_all(i,time,nn) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
        tv_penalty(i,nn) = sum(abs(DiffPhiX_1D(X_hat)),'all');
        awmv_rmse(i,nn) = norm(awmv_all(i,:,nn)-theta_stds1')/norm(theta_stds1);
    end
end
[~,s_i] = min(awmv_rmse);

%% Use Lambda_t and Gamma parameters to do MC
trials = 100;
% Input dirs
dset_name = ['anomaly_noise_MC_Poisson'];

% Output dirs
indep_name = '_indep_ISM';
output_name = '_coupled';
output_subdir = [dset_name,output_name];
indep_subdir = [dset_name,indep_name];

% Setup directories
dataset = fullfile(top_dir,dset_name);
output_dir = fullfile(top_dir,output_subdir);
indep_dir = fullfile(top_dir,indep_subdir);
mkdir(output_dir)

P.params.maxIter = 100;
P.params.rho1 = 1.5;
P.params.rho2 = 0.5;
X_coupled = zeros(N,K,trials,T);
for nn = 1:NN
    indep_data = load(fullfile(indep_dir,[dset_name,'_',num2str(nn),'_','all']));
    X_indep = indep_data.X_indep;
    B = indep_data.B;
    P.set = nn;
    for i = 1:trials 
        P.params.lambda1 = lambda_select(nn,:);
        P.params.lambda2 = gamma_values(nn,s_i(nn));
        X_init = squeeze(X_indep(:,:,i,:));
        [X_hat,~,~,~,~] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,...
                          B(:,:,i),X_init,P.params);
        X_coupled(:,:,i,:) = X_hat;
    end
    save(fullfile(output_dir,[dset_name,'_',num2str(nn),'_','CGTV1']),...
                  'B','X_coupled','P');
end


