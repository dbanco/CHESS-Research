%% MSE for AWMV and characterize noise in AWMV
top_dir = 'D:\CHESS_data\';
T = 50;
M = 20;
noise_std = [0:0.03:0.30];
theta_stds1 = linspace(1,15,T)';
lambda_values = logspace(-4,1,M);
lambda_inds = [15,17,19,21,...
               22,22,22,22,...
               23,23,23];
           
NN = numel(noise_std);
mse = zeros(NN,1);
for ii = 1:NN           
    dset_name = ['singlePeak_noise',num2str(ii)];
    % Output dirs
    output_name = '_indep_ISM';
    output_subdir = [dset_name,output_name];
    output_dir  = fullfile(top_dir,output_subdir);
    i = lambda_inds(ii);
    load(fullfile(output_dir,[dset_name,'_',num2str(i),'_','all']))
    diff = awmv_az - theta_stds1;
    mse(ii) = norm(diff);
end

figure(1)
plot(noise_std,mse,'-o')
ylabel('AWMV MSE')
xlabel('Additive White Gaussian Noise \sigma')