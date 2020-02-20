%% Parameter selection
clear all
close all

% dataset = '/cluster/home/dbanco02/mmpad_polar/ring1_zero/';
% output_dir = '/cluster/shared/dbanco02/mmpad_1D_indep_param_1/';

% output_dir = 'D:\CHESS_data\mmpad_1D_coupled_simul_';
output_dir = 'D:\CHESS_data\noise_1D_coupled_simul_4_';
num_ims = 10;
baseFileName = 'fista_fit_%i_%i.mat';

% Gamma values
gamma_vals = [0.00025,0.0005,0.00075,0.001,0.0025,0.005,0.0075,...
              0.01,0.02,0.03 0.04 0.05,0.08,0.1,0.15,0.2,0.25,...
              0.3,0.4,0.5,0.75,1,2,5,10]; 
M = numel(gamma_vals);

% Universal Parameters
% Ring sampling parameters
load([output_dir,'1a\',sprintf(baseFileName,1,1)])

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
end
A0 = unshifted_basis_vector_stack_norm2(P);

% Construct distance matrix
Threshold = 32;
D = constructDistanceMatrix_1D(P,Threshold);

% Get error, sparsity, awmv
err_select = zeros(M,num_ims);
l0_select = zeros(M,num_ims);
l1_select = zeros(M,num_ims);
awmv_az = zeros(M,num_ims);
vdfs = zeros(P.num_var_t,M,num_ims);
obj1 = zeros(M,1);
obj2 = zeros(M,1);
obj3 = zeros(M,1);

for k = 1:M
    fprintf('%i of %i\n',k,M)
    for j = 1:num_ims
        load(fullfile([output_dir,num2str(k),'aa\'],sprintf(baseFileName,1,j)))
        
        % Fit objective
        fit = Ax_ft_1D(A0ft_stack,x_hat);
        b = polar_image./norm(polar_image(:));
        
        err_select(k,j) = norm(b(:)-fit(:));
        obj1(k) = obj1(k) + 0.5*norm(b(:)-fit(:))^2;
        
        % Sparsity objective
        l0_select(k,j) = sum(x_hat(:)>0);
        l1_select(k,j) = sum(x_hat(:));
        obj2(k) = obj2(k) + l1_select(k,j)*P.params.lambda;
       
        % Vdf
        az_signal = squeeze(sum(x_hat,1));
        var_sum = sum(az_signal(:));
        vdfs(:,k,j) = az_signal./var_sum;
        awmv_az(k,j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end
end

% Compute Wasserstein objective
%Pc.wLam
for k = 1:M
    fprintf('Wasserstein %i of %i\n',k,M)
    for j = 1:num_ims-1
        wass_dist = WassersteinObjective(vdfs(:,k,j),{vdfs(:,k,j+1)},Pc.wLam,D);
        obj3(k) = obj3(k) + wass_dist;
    end
end

% Load statistics for independently fit data
init_dir = 'D:\CHESS_data\noise_1D_coupled_simul_4_init4aa';
awmv_az_init = zeros(num_ims,1);
err_indep = zeros(num_ims,1);
l0_indep = zeros(num_ims,1);
l1_indep = zeros(num_ims,1);
vdfs_indep = zeros(P.num_var_t,num_ims);

for j = 1:num_ims
    load(fullfile(init_dir,sprintf(baseFileName,1,j)))
    
    % Fit objective
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    b = polar_image./norm(polar_image(:));
    err_indep(j) = norm(b(:)-fit(:));
    
    % Sparsity objective
    l0_indep(j) = sum(x_hat(:)>0);
    l1_indep(j) = sum(x_hat(:));

    % Vdf
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    vdfs_indep(:,j) = az_signal./var_sum;
    awmv_az_init(j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
end

wass_indep = 0;
for j = 1:num_ims-1
    wass_dist = WassersteinObjective(vdfs_indep(:,j),{vdfs_indep(:,j+1)},Pc.wLam,D);
    wass_indep = wass_indep + wass_dist;
end

%% Plot error, std, difference from indep awmv
close all
% Plot AWMV
figure(1)
legend_str = cell(M+1,1);
legend_str{1} = 'indep';
hold on
plot(awmv_az_init,'-')
for k = 1:M
    hold on
    plot(awmv_az(k,:),'-')
    legend_str{k+1} = sprintf('%0.04f',gamma_vals(k));
end

ylabel('AWMV')
xlabel('t')
legend(legend_str,'location','best')

[sort_gamma, sort_i] = sort(gamma_vals);

% Plot err
mean_err = mean(err_select,2);
figure(2)
hold on
plot(0,mean(err_indep),'o')
plot(sort_gamma,mean_err(sort_i),'o-')
ylabel('Average Error')
xlabel('Coupling parameter')

% Plot l1 norm
mean_l1 = mean(l1_select,2);
figure(3)
hold on
plot(0,mean(l1_indep),'o')
plot(sort_gamma,mean_l1(sort_i),'o-')
ylabel('l1-norm')
xlabel('Coupling parameter')

% Plot number nonzeros coefficients
mean_l0 = mean(l0_select,2);
figure(4)
hold on
plot(0,mean(l0_indep),'o')
plot(sort_gamma,mean_l0(sort_i),'o-')
ylabel('l0-norm')
xlabel('Coupling parameter')

% Plot wass dist
wass_total = obj3;
figure(44)
hold on
plot(0,wass_indep,'o')
plot(sort_gamma,wass_total(sort_i),'o-')
ylabel('Wasserstein distance')
xlabel('Coupling parameter')

% Plot AWMV
figure(5)
legend_str = {};
legend_str{1} = 'indep';
hold on
plot(awmv_az_init,'-')
kk = 2;
for k = [4,5,M]
    hold on
    plot(awmv_az(sort_i(k),:),'-')
    legend_str{kk} = sprintf('%0.04f',sort_gamma(k));
    kk = kk + 1;
end
legend(legend_str,'Location','Best')

% Plot L-curve 1
figure(6)
plot(obj2(sort_i),obj1(sort_i),'o-')
xlabel('l1-norm')
ylabel('Error')

% Plot L-curve 2
figure(7)
plot(obj3(sort_i),obj1(sort_i),'o-')
xlabel('Wasserstein distance')
ylabel('Error')

% Plot L-curve 3
figure(8)
plot(obj3(sort_i),obj2(sort_i),'o-')
xlabel('Wasserstein distance')
ylabel('l1-norm')

% Plot L-curve 4
figure(9)
plot(obj3(sort_i),obj1(sort_i)+obj2(sort_i),'o-')
xlabel('Wasserstein distance')
ylabel('Error+l1-norm')

%% Plot fits

figure(222)
[ha2, pos2] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az = zeros(num_ims,1);
im_ind = 1;
trial_k = 4;
for image_num = 1:10
    
    load(fullfile([output_dir,num2str(sort_i(trial_k)),'a\'],sprintf(baseFileName,1,image_num)))
    
    polar_vector = polar_image
    fit = Ax_ft_1D(A0ft_stack,x_hat)*norm(polar_vector(:));
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = zeroPad(polar_vector,P.params.zeroPad);
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x_hat(:)>1e-6)),'location','northeast')
    im_ind = im_ind + 1;
end

%% Plot vdfs
figure(333)
[ha3, pos3] = tight_subplot(5,10,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az = zeros(num_ims,1);
im_ind = 1;
trial_k = 4;
for image_num = 1:5:250
    
    load(fullfile([output_dir,num2str(sort_i(trial_k)),'a\'],sprintf(baseFileName,1,image_num)))
    
    polar_vector = squeeze(sum(polar_image,1));
    fit = Ax_ft_1D(A0ft_stack,x_hat)*norm(polar_vector(:));
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    % Plot
    axes(ha3(im_ind))
    plot(az_signal/var_sum)
    legend(sprintf('%i',sum(x_hat(:)>1e-6)),'location','northeast')
    im_ind = im_ind + 1;
end

%% Plot vdfs
figure(333)
[ha3, pos3] = tight_subplot(5,10,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az = zeros(num_ims,1);
im_ind = 1;
trial_k = 4;
for image_num = 1:5:250
    
    load(fullfile([output_dir,'init\'],sprintf(baseFileName,1,image_num)))
    
    polar_vector = squeeze(sum(polar_image,1));
    fit = Ax_ft_1D(A0ft_stack,x_hat)*norm(polar_vector(:));
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    % Plot
    axes(ha3(im_ind))
    plot(az_signal/var_sum)
    legend(sprintf('%i',sum(x_hat(:)>1e-6)),'location','northeast')
    im_ind = im_ind + 1;
end

%% Compute selection criteria
lambda_indices = zeros(num_ims,1);
noise_est = zeros(num_ims,1);
norm_data = zeros(num_ims,1);
noise_eta = zeros(num_ims,1);
norm_ratio = zeros(num_ims,1);
for image_num = 1:num_ims
    fprintf('%i of %i\n', image_num, num_ims)
    im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
    % Zero pad image
    b = squeeze(sum(im_data.polar_image,1));
    % Scale image by 2-norm
    bn = b/norm(b(:));
    
    [~,cN] = size(b);
    kernel = shift1D(A0(:,1),round(cN/2));
%     kernel = kernel(1:35,:);
    bn_hat = conv(bn,kernel,'same');
    noise_est(image_num) = norm(bn-bn_hat);
    norm_data(image_num) = norm(b);
    norm_ratio(image_num) = norm(b)/norm(b,1);
end

%%  Select paramters 

% norm scaling method
% noise_eta = norm_data./max(norm_data).*max(noise_est);

% purely residual based method
noise_eta = 0.1;

% some other method
% noise_eta = max(noise_est)*norm_ratio;

% Criterion
discrep_crit = abs(err_select'-noise_eta);
% discrep_crit = abs(l0_select'-80);
% discrep_crit = abs(err_select'-repmat(noise_eta,1,N));



[lambda_indices,~] = find(discrep_crit' == min(discrep_crit'));
param_select = lambda_vals(lambda_indices);

figure(4)
plot(noise_eta)
title('Noise estimate')

figure(3)
plot(noise_est)
title('norm(b-b_t)')

figure(2)
plot(param_select)
title('Paramters selected')

figure(6)
plot_err_select = err_select;
plot_err_select(err_select > 100) = 1;

% plot(plot_err_select)
loglog(repmat(lambda_vals',[1,500]),plot_err_select)
hold on
for i = 1:500
   plot(param_select(i),err_select(lambda_indices(i),i),'o') 
end
title('Error')
%% Load with selected parameters and plot

figure(222)
[ha2, pos2] = tight_subplot(5,10,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az = zeros(num_ims,1);
im_ind = 1;
for image_num = 1:5:250
    
    load(fullfile(output_dir,sprintf(baseFileName,lambda_indices(image_num),image_num)))
    
    polar_vector = squeeze(sum(polar_image,1));
    fit = Ax_ft_1D(A0ft_stack,x_hat)*norm(polar_vector(:));
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    
    
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(polar_vector)
    plot(fit(51:end))
    legend(sprintf('%i',sum(x_hat(:)>1e-6)),'location','northeast')
    im_ind = im_ind + 1;
end

%% Plot resulting AWMV
for image_num = 1:500
    load(fullfile(output_dir,sprintf(baseFileName,lambda_indices(image_num),image_num)))
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
end

figure(1)
hold on
plot(awmv_az,'-')
legend('Fit Discrep','Truth','Fit Single \lambda','Location','best')
ylabel('AWMV')
xlabel('t')

%% Plot basis functions
figure(5)
for i = 1:P.num_var_t
   kernel = shift1D(A0(:,i),round(cN/2)+50);
   hold on
   plot(kernel)
end

