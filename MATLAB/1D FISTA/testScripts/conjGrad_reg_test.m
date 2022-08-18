%% Parameter selection
clear all
close all
P.set = 1;

dataset = 'D:\CHESS_data\simulated_two_spot_1D_gnoise4_nonorm_3';

num_ims = 20;
dset_name = 'gnoise4_nonorm';

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta = size(polar_vector,1);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 20;
P.var_theta = linspace(P.dtheta/2,50,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% lasso/admm params
params.lambda1 = 0.0036;
params.lambda2 = 0.0001;
params.rho1 = 1;
params.rho2 = 1;

params.tau = 1.05;
params.mu = 2;
params.adaptRho = 1;
params.alpha = 1.8;
params.maxIter = 800;
params.conjGradIter = 50;

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.tolerance = 1e-6;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.plotProgress = 0;
params.verbose = 1;

P.params = params;
   
baseFileName = 'fista_fit_%i_%i.mat';

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
A0 = unshifted_basis_vector_stack_norm2_zpad(P);

% Run grid search ISM
T = 20;
[N,M] = size(A0ft_stack);
X_init = zeros(N,M,T);
B = zeros(N,T);
for image_num = 1:T
    im_data = load(fullfile(dataset,[prefix,'_',num2str((image_num)),'.mat']));
    % Zero pad image
    B(:,image_num) = im_data.polar_vector;
end

[X_hat_tvx, err, obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_TVx_1D(A0ft_stack,B,X_init,params);
% [X_hat1, err, obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_1D(A0ft_stack,B,X_init,params);
%% Plot results
% for t = 1:T
%     x = X_hat(:,:,t);
%     x(x < sum(x(:),'all')*1e-4) = 0;
%     bb = Ax_ft_1D(A0ft_stack,x);
%     figure(1)
%     subplot(T/2,2,t)
%     plot(B(:,t))
%     hold on
%     plot(bb)  
%     l0_norm = sum((x > 0),'all');
%     legend('Data',sprintf('Fit: %i',l0_norm))
%     
%     figure(3)
%     subplot(T/2,2,t)
%     imagesc(x)
% end
% 
% for t = 1:(T-1)
%     temp = sum( ( X_hat(:,:,t)-X_hat(:,:,t+1) ).^2 );
% end
% 
% 
% sum(abs(DiffX_1D(X_hat)),'all')
% figure(2)
% subplot(2,2,1)
% plot(obj(2:end))
% title('Objective')
% 
% subplot(2,2,2)
% plot(err(2:end))
% title('Error')
% 
% subplot(2,2,3)
% plot(l1_norm(2:end))
% title('l1 norm')
% 
% subplot(2,2,4)
% plot(tv_penalty(2:end))
% title('TVx')


%% Plot vdf and compare to independent fit
n_level = 3;
dataset_num = num2str(n_level);
num_ims = 20;

top_dir = 'D:\CHESS_data';

% Input dirs
dset_name = 'gnoise4_nonorm';
dset_subdir = ['simulated_two_spot_1D_',dset_name,'_',num2str(n_level)];
indep_name = 'ADMM_CG_indep1';
indep_subdir = [dset_subdir,'_indep'];
init_subdir =  [dset_subdir,'_simul_init'];

% Output dirs
output_name = 'gnoise4_nonorm_coupled_CG_TVx3';
output_subdir = [dset_subdir,'_coupled'];

% Setup directories
dataset =  fullfile(top_dir,dset_subdir);
indep_dir = fullfile(top_dir,indep_name,indep_subdir);
init_dir =  fullfile(top_dir,indep_name,init_subdir);
output_dir  = fullfile(output_name,output_subdir);
baseFileName = 'fista_fit_%i.mat';

ind_data = load(fullfile(indep_dir,sprintf(baseFileName,10)));
% ind_data = load(fullfile(init_dir,sprintf(baseFileName,1)));

ind_data.P.params

figure(4)
vdf_tvx = squeeze(sum(X_hat_tvx,1))';
% vdf_indep1 = squeeze(sum(X_hat1,1))';
vdf_indep2 = squeeze(sum(ind_data.X_hat,1))';

subplot(3,1,1)
imagesc(vdf_tvx)
title('TVx, lambda2=0.0001')

% subplot(3,1,2)
% imagesc(vdf_indep1)
% title('Indep1, lambda1=.0036')

subplot(3,1,3)
imagesc(vdf_indep2)
title('Indep2, lambda1=0.0036')

