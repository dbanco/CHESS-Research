%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\CHESS_data\';
% top_dir = '/cluster/shared/dbanco02';\

noise_std = [0:0.03:0.30];
% lambda_inds = [15,17,19,21,...
%                22,22,22,23,...
%                23,23,24];
lambda_inds = [15,17,19,21,...
               22,22,22,23,...
               23,23,24];
NN = numel(noise_std);
ii = 9;
% Input dirs
dset_name = ['singlePeak_noise',num2str(ii)];

% Output dirs
output_name = '_coupled_CGTV';
output_subdir = [dset_name,output_name];


% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
mkdir(output_dir)

num_ims = 50;

% File Parameters
P.baseFileName = 'couple_fit_%i_%i.mat';
P.dataset = dataset;

% Data/Dictionary Parameters
% Zero padding and mask


N = 101;
K = 20;
M = 30;
T = num_ims;
zPad = 0;
zMask = [];

P.dataScale = 1;
P.lambda_values = logspace(-4,1,M);
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
P.params.maxIter = 100;
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

% Load data
theta_stds1 = linspace(1,15,T);
B = zeros(N+2*zPad,T);
for j = 1:T
    b = gaussian_basis_1D( N, N/2, theta_stds1(j)^2) + randn(N,1)*noise_std(ii);

    B(:,j) = b;
end

%% Independent Solution
x_init = zeros(N,K);
X_indep = zeros(N,K,T);
rho_out = zeros(T,1);
for i = lambda_inds(ii)
    P.set = i;
    P.params.lambda1 = P.lambda_values(i);
    for t = 1:T
        % Solve
        [x_hat,obj,err,l1_norm,rho] = convADMM_LASSO_Sherman_1D(A0ft_stack,B(:,t),x_init,P.params);  
        X_indep(:,:,t) = x_hat;
        rho_out(t) = rho;
    end
end

%% CG Indep
X_indep_CG = zeros(N,K,T);
P.params.maxIter = 200;
for i = lambda_inds(ii)
    P.set = i;
    P.params.lambda1 = P.lambda_values(i);
    for t = 1:T
        % Solve
        [x_hat,obj,err,l1_norm,~] = convADMM_LASSO_CG_1D(A0ft_stack,B(:,t),x_init,P.params);  
        X_indep_CG(:,:,t) = x_hat;

    end
end

%% Coupled Solution
P.params.rho1 = 1;%min(rho_out);
P.params.rho2 = 1;
P.params.lambda2 = 0.1;

for i = lambda_inds(ii)
    P.set = i;
    P.params.lambda1 = ones(T,1)*P.lambda_values(i);
    [X_hat,err,obj,l1_norm,tv_penalty] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,B,zeros(N,K,T),P.params);
    save(fullfile(output_dir,sprintf(P.baseFileName,P.set,t)),'X_indep','X_hat','err','obj','l1_norm','tv_penalty','P');
end
%% Plot fit 
close all
fits_fig = figure(1);
[ha2, ~] = tight_subplot(10,5,[.005 .005],[.01 .01],[.01 .01]);
fits_fig_indep = figure(11);
 [ha11, ~] = tight_subplot(10,5,[.005 .005],[.01 .01],[.01 .01]);
fits_fig_indep_CG = figure(111);
 [ha111, ~] = tight_subplot(10,5,[.005 .005],[.01 .01],[.01 .01]);
 
awmv_az_indep = zeros(T,1);
awmv_az_indep_CG = zeros(T,1);
awmv_az = zeros(T,1);
rse_coup = zeros(T,1);
rse_indep = zeros(T,1);
rse_indep_CG = zeros(T,1);
for t = 1:T
    x = X_hat(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x);
    rse_coup(t) = norm(fit-B(:,t))/norm(B(:,t));
    az_signal = squeeze(sum(x,1));
    var_sum = squeeze(sum(az_signal(:)));
    awmv_az(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    axes(ha2(t))
    hold on
    plot(B(:,t))
    plot(fit)
    
    x = X_indep(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x);
    rse_indep(t) = norm(fit-B(:,t))/norm(B(:,t));
    az_signal = squeeze(sum(x,1));
    var_sum = squeeze(sum(az_signal(:)));
    awmv_az_indep(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    axes(ha11(t))
    hold on
    plot(B(:,t))
    plot(fit)
    
    x = X_indep_CG(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x);
    rse_indep_CG(t) = norm(fit-B(:,t))/norm(B(:,t));
    az_signal = squeeze(sum(x,1));
    var_sum = squeeze(sum(az_signal(:)));
    awmv_az_indep_CG(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    axes(ha111(t))
    hold on
    plot(B(:,t))
    plot(fit)

end
legend('fit','b')
save(fullfile(output_dir,[dset_name,'_',num2str(P.set),'_','all']),...
        'B','X_hat','X_indep','X_indep_CG',...
        'awmv_az','P');
    
%% Plot awmv, rel err 
awmv_fig = figure(2);
hold on
plot(theta_stds1,theta_stds1,'-','Linewidth',2)
% plot(theta_stds1,awmv_az_indep,'Linewidth',1)
plot(theta_stds1,awmv_az_indep_CG,'Linewidth',1)
plot(theta_stds1,awmv_az,'Linewidth',1)
xlabel('\sigma')
ylabel('AWMV')

mse1 = norm(theta_stds1-awmv_az_indep);
mse2 = norm(theta_stds1-awmv_az_indep_CG);
mse3 = norm(theta_stds1-awmv_az);

% mse1 = norm(theta_stds1-awmv_az_indep)/norm(theta_stds1);
% mse2 = norm(theta_stds1-awmv_az_indep_CG)/norm(theta_stds1);
% mse3 = norm(theta_stds1-awmv_az)/norm(theta_stds1);

legend('truth',['indep- ',num2str(mse2)],...
                ['coupled- ',num2str(mse3)],...
                'location','best')
            
% legend('truth',['indep- ',num2str(mse1)],...
%                 ['indep_{CG}- ',num2str(mse2)],...
%                 ['coupled- ',num2str(mse3)],...
%                 'location','best')

rel_err_fig = figure(3);
hold on
plot(theta_stds1,rse_indep,'-','Linewidth',2)
plot(theta_stds1,rse_indep_CG,'-','Linewidth',2)
plot(theta_stds1,rse_coup,'-','Linewidth',1)
legend('indep','indep_{CG}','coupled','location','best')
%end

%% Check awmv
% close all
% ii = 5;
% i = lambda_inds(ii);
% dset_name = ['singlePeak_noise',num2str(ii)];
% % Output dirs
% output_name = '_indep_ISM';
% output_subdir = [dset_name,output_name];
% output_dir  = fullfile(top_dir,output_subdir);
% 
% load(fullfile(output_dir,[dset_name,'_',num2str(i),'_','all']))
% awmv_fig = figure(2);
% hold on
% plot(theta_stds1,theta_stds1,'-','Linewidth',2)
% plot(theta_stds1,awmv_az,'Linewidth',1)
% xlabel('\sigma')
% ylabel('AWMV')
% awmv_fig.Position = [800 700 200 200];
%     set(gca,'XTickLabel',[]);
%     set(gca,'YTickLabel',[]);
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(gca,'XColor', 'none','YColor','none')

%% Produce 3 example fits
% fits_fig = figure(1);
% [ha2, ~] = tight_subplot(1,3,[.005 .005],[.01 .01],[.01 .01]); 
% lambda_inds = [15,17,19,21,...
%                22,22,22,23,...
%                23,23,24];
% ii = 3;
% i = lambda_inds(ii);
% dset_name = ['singlePeak_noise',num2str(ii)];
% % Output dirs
% output_name = '_indep_ISM';
% output_subdir = [dset_name,output_name];
% output_dir  = fullfile(top_dir,output_subdir);
% awmv_az = zeros(T,1);
% load(fullfile(output_dir,[dset_name,'_',num2str(i),'_','all']))
% A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
% im_ind=1;
% for t = [1,25,50]
%     x = X_hat(:,:,t);
%     fit = Ax_ft_1D(A0ft_stack,x);
%     az_signal = squeeze(sum(x,1));
%     var_sum = squeeze(sum(az_signal(:)));
%     awmv_az(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
%     
%     axes(ha2(im_ind))
%     hold on
%     plot(B(:,t),'Linewidth',2)
%     plot(fit,'Linewidth',1)
%     im_ind = im_ind+1;
%     set(gca,'XTickLabel',[]);
%     set(gca,'YTickLabel',[]);
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(gca,'XColor', 'none','YColor','none')
% end
% fits_fig.Position = [800 800 300 100];


