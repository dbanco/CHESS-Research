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
output_name = '_indep_ISM';
output_subdir = [dset_name,output_name];
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
load(fullfile(output_dir,[dset_name,'_',num2str(1),'_','all']))
[N,K,M,T] = size(X_indep);
theta_stds1 = [7*ones(1,T/2),12*ones(1,T/2)];

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Analysis of MC runs
NN = 11;
awmv_indep = zeros(NN,M,T);
rse_indep = zeros(NN,M,T);
mse_indep = zeros(NN,M,T);
l1_norm = zeros(NN,M,T);
% x_avg = zeros(NN,N,K,T);
fit_avg = zeros(NN,N,T);
fit_1 = zeros(NN,N,T);
awmv_rse = zeros(NN,M);

% Noise levels
for nn = 1:11
    output_subdir = [dset_name,output_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    load(fullfile(output_dir,[dset_name,'_',num2str(nn),'_','all']))
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
            l1_norm(nn,i,time) = sum(abs(x(:)));
            mse_indep(nn,i,time) = norm(fit-B(:,time));
            rse_indep(nn,i,time) = norm(fit-B(:,time))/norm(B(:,time));
%             x_avg(:,:,time) = x_avg(:,:,time) + squeeze(x)/100;
            fit_avg(nn,:,time) = squeeze(fit_avg(nn,:,time))' + squeeze(fit)/100;
        end
        awmv_rse(nn,i) = norm( squeeze(awmv_indep(nn,i,:))'-...
                         theta_stds1 )/norm(theta_stds1);
    end
end

%% Waterfall for highest noise level
figure(1)
nn = 11;
fit_nn_1 = squeeze(fit_1(11,:,:));
fit_nn = squeeze(fit_avg(11,:,:));

subplot(1,3,1)
waterfall(fit_nn')
title('Averaged Fits for 100 trials')

subplot(1,3,2)
waterfall(fit_nn_1')
title('Single Fit')

subplot(1,3,3)
waterfall(B')
title('Data')

%% Fit RMSE per NN
figure(2)
avg_rse_indep = squeeze(mean(rse_indep,2));
rse_nn = mean(avg_rse_indep,2);
plot(noise_std,rse_nn)
xlabel('v')
title('Fit Error')
ylabel('Average Fit RSE over 100 trials')

%% AWMV RMSE per NN
figure(3)
avg_awmv = squeeze(mean(awmv_indep,2));
avg_awmv_rse = zeros(NN,1);
mean_awmv_rse = mean(awmv_rse,2);
for nn = 1:NN
    avg_awmv_rse(nn) = norm(avg_awmv(nn,:)-theta_stds1)/norm(theta_stds1);
end
plot(noise_std,mean_awmv_rse)
title('AWMV Error')
xlabel('v')
ylabel('Average AWMV RSE over 100 trials')

figure(4)
plot(avg_awmv(10,:)')
