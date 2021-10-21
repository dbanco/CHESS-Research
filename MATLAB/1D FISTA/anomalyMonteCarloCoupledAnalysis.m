%% Analysis of MC runs
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\CHESS_data\';
% top_dir = '/cluster/shared/dbanco02';

noise_std = [0:0.03:0.30];
NN = numel(noise_std);

% Load in Parameters
dset_name = ['anomaly_noise_MC'];
output_name = '_coupled';
output_subdir = [dset_name,output_name];
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
load(fullfile(output_dir,[dset_name,'_',num2str(1),'_','CGTV1']))
[N,K,M,T] = size(X_coupled);
theta_stds1 = [7*ones(1,T/2),12*ones(1,T/2)]';

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Analysis of MC runs
NN = 11;
awmv_coup = zeros(NN,M,T);
rse_coup = zeros(NN,M,T);
mse_coup = zeros(NN,M,T);
l1_norm_coup = zeros(NN,M,T);
% x_avg = zeros(NN,N,K,T);
fit_avg_coup = zeros(NN,N,T);
fit_1 = zeros(NN,N,T);
awmv_rse_coup = zeros(NN,M);

%% Noise levels
for nn = 1:11
    output_subdir = [dset_name,output_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    load(fullfile(output_dir,[dset_name,'_',num2str(nn),'_','CGTV1']))
    % 100 trials
    for i = 1:M
        for time = 1:T
            x = X_coupled(:,:,i,time);
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_coup(nn,i,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
            fit = Ax_ft_1D(A0ft_stack,x);
            if( i == 2)
                fit_1(nn,:,time) = fit;
            end
            l1_norm_coup(nn,i,time) = sum(abs(x(:)));
            mse_coup(nn,i,time) = norm(fit-B(:,time,i));
            rse_coup(nn,i,time) = norm(fit-B(:,time,i))/norm(B(:,time,i));
            fit_avg_coup(nn,:,time) = squeeze(fit_avg_coup(nn,:,time))' + squeeze(fit)/100;
        end
        awmv_rse_coup(nn,i) = norm( squeeze(awmv_coup(nn,i,:))-...
                         theta_stds1 )/norm(theta_stds1);
    end
end

%% Analyze independent solutions
indep_name = '_indep_ISM';
indep_subdir = [dset_name,indep_name];
indep_dir  = fullfile(top_dir,indep_subdir);
load(fullfile(indep_dir,[dset_name,'_',num2str(1),'_','all']))

awmv_indep = zeros(NN,M,T);
rse_indep = zeros(NN,M,T);
mse_indep = zeros(NN,M,T);
l1_norm_indep = zeros(NN,M,T);
% x_avg = zeros(NN,N,K,T);
fit_avg_indep = zeros(NN,N,T);
awmv_rse_indep = zeros(NN,M);

for nn = 1:11
    load(fullfile(indep_dir,[dset_name,'_',num2str(nn),'_','all']))
    % 100 trials
    for i = 1:M
        for time = 1:T
            x = X_indep(:,:,i,time);
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_indep(nn,i,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
            fit = Ax_ft_1D(A0ft_stack,x);
            if( i == 2)
                fit_1(nn,:,time) = fit;
            end
            l1_norm_indep(nn,i,time) = sum(abs(x(:)));
            mse_indep(nn,i,time) = norm(fit-B(:,time,i));
            rse_indep(nn,i,time) = norm(fit-B(:,time,i))/norm(B(:,time,i));
            fit_avg_indep(nn,:,time) = squeeze(fit_avg_indep(nn,:,time))' + squeeze(fit)/100;
        end
        awmv_rse_indep(nn,i) = norm( squeeze(awmv_indep(nn,i,:))-...
                         theta_stds1 )/norm(theta_stds1);
    end
end

%% Inspect AWMV curves

figure(21)
hold on
nn = 8;
% trial = 3;
% errorbar(squeeze(mean(awmv_coup(nn,:,:),2)),...
%          squeeze(std(awmv_coup(nn,:,:),1,2)) )
errorbar(squeeze(mean(awmv_indep(nn,:,:),2)),...
         squeeze(std(awmv_indep(nn,:,:),1,2)) )
plot(theta_stds1)
% legend('Coupled','Indep','Truth','Location','Best')


figure(22)
hold on
nn = 8;
% trial = 3;
errorbar(squeeze(mean(awmv_coup(nn,:,:),2)),...
         squeeze(std(awmv_coup(nn,:,:),1,2)) )
% errorbar(squeeze(mean(awmv_indep(nn,:,:),2)),...
%          squeeze(std(awmv_indep(nn,:,:),1,2)) )
plot(theta_stds1)
% legend('Coupled','Indep','Truth','Location','Best')

figure(23)
hold on
nn = 8;
trial = 15;
plot(squeeze(awmv_coup(nn,trial,:)))
plot(squeeze(awmv_indep(nn,trial,:)))
plot(theta_stds1)
legend('Coupled','Indep','Truth','Location','Best')


%% Waterfall for highest noise level compare
figure(1)
nn = 11;
fit_nn_1 = squeeze(fit_1(11,:,:));
fit_nn = squeeze(fit_avg_coup(11,:,:));

subplot(1,3,1)
waterfall(fit_nn')
title('Averaged Fits for 100 trials')

subplot(1,3,2)
waterfall(fit_nn_1')
title('Single Fit')

subplot(1,3,3)
waterfall(B(:,:,1)')
title('Data')

%% Fit RMSE per NN
figure(2)
avg_rse_indep = squeeze(mean(rse_coup,2));
rse_nn = mean(avg_rse_indep,2);
plot(noise_std,rse_nn)
xlabel('v')
title('Fit Error')
ylabel('Average Fit RSE over 100 trials')

%% AWMV RMSE per NN comparison
f3 = figure(3);
avg_awmv_coup =  squeeze(mean(awmv_rse_coup,2));
avg_awmv_indep = squeeze(mean(awmv_rse_indep,2));

std_awmv_coup =  squeeze(std(awmv_rse_coup,1,2));
std_awmv_indep = squeeze(std(awmv_rse_indep,1,2));

hold on
errorbar(noise_std,avg_awmv_indep,std_awmv_indep,'LineWidth',2,'color','blue')
errorbar(noise_std,avg_awmv_coup,std_awmv_coup,'LineWidth',2,'color','red')
xlabel('v (noise variance)','FontSize',16)
ylabel('Average AWMV RSE (100 trials)','FontSize',16)
legend('\gamma = 0','\gamma = \gamma*','Location','NorthWest','FontSize',16)
grid on
f3.Position = [1000,500,700 350];

% figure(4)
% plot(avg_awmv(10,:)')
