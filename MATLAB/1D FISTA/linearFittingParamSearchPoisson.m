%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\CHESS_data\';
% top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = ['singlePeak_noise_Poisson'];

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

% Load data
theta_stds1 = linspace(1,15,T);
B = zeros(N,T);
for t = 1:T
    close all
    b = gaussian_basis_1D( N, N/2, theta_stds1(t)^2);
    rms = sqrt(sum(b.^2)/N);
    b = b*100; % + randn(N,1)*noise_std(nn);
    bn = poissrnd(b);
    B(:,t) = bn;
    plot(b);hold on;plot(bn);
    
end

%% Independent Solution
% x_init = zeros(N,K);
% X_indep = zeros(N,K,M,T);
% for i = 1:M
%     P.set = i;
%     P.params.lambda1 = P.lambda_values(i);
%     for t = 1:T
%         % Solve
%         [x_hat,obj,err,l1_norm,~] = convADMM_LASSO_Sherman_1D(A0ft_stack,B(:,t),x_init,P.params);  
%         X_indep(:,:,i,t) = x_hat;
%     end
% end
% save(fullfile(output_dir,[dset_name,'_',num2str(P.set),'_','all']),...
%         'B','X_indep','P');


%% Plot fit
dset_name = ['singlePeak_noise_Poisson'];
output_name = '_indep_ISM';
output_subdir = [dset_name,output_name];
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
load(fullfile(output_dir,[dset_name,'_',num2str(50),'_','all']))

close all
fits_fig_indep = figure(11);
 [ha11, ~] = tight_subplot(10,5,[.005 .005],[.01 .01],[.01 .01]);
awmv_az_indep = zeros(T,1);
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
%     axes(ha11(time))
%     hold on
%     plot(B(:,time))
%     plot(fit)
    end
end

% Plot selected awmv
close all
select_ind = zeros(T,1);
for time = 1:T
    crit = abs(l1_norm(:,time)*0.6).^2 + abs(mse_indep(:,time)).^2;
    select_ind(time) = find( (crit == min(crit)),1 );
    x = X_indep(:,:,select_ind(time),time);
    az_signal = squeeze(sum(x,1));
    var_sum = squeeze(sum(az_signal(:)));
    awmv_az_indep(time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
end
select_ind
figure(1)
hold on
plot(awmv_az_indep)
plot(theta_stds1)
legend('awmv','truth')
%% Plot L-curve
figure(3)
plot(l1_norm,mse_indep,'o-')
% waterfall(1:T,l1_norm,rse_indep)
xlabel('l_1 norm')
ylabel('Rel Error')

%% Plot L-curve
nn = 1
figure(4)
plot(l1_norm(15:50,nn),mse_indep(15:50,nn),'o-')
% waterfall(1:T,l1_norm,rse_indep)
xlabel('l_1 norm')
ylabel('Rel Error')


%% Plot awmv, rel err 
% awmv_fig = figure(2);
% hold on
% plot(theta_stds1,theta_stds1,'-','Linewidth',2)
% % plot(theta_stds1,awmv_az_indep,'Linewidth',1)
% plot(theta_stds1,awmv_az_indep_CG,'Linewidth',1)
% plot(theta_stds1,awmv_az,'Linewidth',1)
% xlabel('\sigma')
% ylabel('AWMV')
% 
% mse1 = norm(theta_stds1-awmv_az_indep);
% mse2 = norm(theta_stds1-awmv_az_indep_CG);
% mse3 = norm(theta_stds1-awmv_az);
% 
% % mse1 = norm(theta_stds1-awmv_az_indep)/norm(theta_stds1);
% % mse2 = norm(theta_stds1-awmv_az_indep_CG)/norm(theta_stds1);
% % mse3 = norm(theta_stds1-awmv_az)/norm(theta_stds1);
% 
% legend('truth',['indep- ',num2str(mse2)],...
%                 ['coupled- ',num2str(mse3)],...
%                 'location','best')
%             
% % legend('truth',['indep- ',num2str(mse1)],...
% %                 ['indep_{CG}- ',num2str(mse2)],...
% %                 ['coupled- ',num2str(mse3)],...
% %                 'location','best')
% 
% rel_err_fig = figure(3);
% hold on
% plot(theta_stds1,rse_indep,'-','Linewidth',2)
% plot(theta_stds1,rse_indep_CG,'-','Linewidth',2)
% plot(theta_stds1,rse_coup,'-','Linewidth',1)
% legend('indep','indep_{CG}','coupled','location','best')
% %end

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

