%% Parameter selection
clear all
close all

% dataset = '/cluster/home/dbanco02/mmpad_polar/ring1_zero/';
% output_dir = '/cluster/shared/dbanco02/mmpad_1D_indep_param_1/';

dataset = 'D:\MMPAD_data\ring1_zero\';
output_dir = 'D:\CHESS_data\mmpad_1D_indep_param_5\';
num_ims = 500;
baseFileName = 'fista_fit_%i_%i.mat';

% Lambda values
lambda_vals = logspace(-3,1,30);
N = numel(lambda_vals);

% Universal Parameters
% Ring sampling parameters
prefix = 'mmpad_img';
load([output_dir,sprintf(baseFileName,1,1)])
polar_image = squeeze(sum(polar_image,1));
polar_image = zeroPad(polar_image,P.params.zeroPad);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
end
A0 = unshifted_basis_vector_stack_norm2(P);

% Get error and sparsity
err_select = zeros(N,num_ims);
l0_select = zeros(N,num_ims);
for i = 1:N
    fprintf('%i of %i\n',i,N)
    for j = 1:num_ims
        load(fullfile(output_dir,sprintf(baseFileName,i,j)))
        err_select(i,j) = err(end);
        l0_select(i,j) = sum(x_hat(:)>0);
    end
end

% function save_output(output_dir,baseFileName,x_hat,err,polar_image,P,image_num)
%     save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','polar_image','P')
% end
% function save_obj(output_dir,pass,image_num,obj)
%     save(fullfile(output_dir,sprintf('objective_%i_%i.mat',pass,image_num)),'obj')
% end

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
   kernel = shift1D(A0(:,i),round(cN/2));
   hold on
   plot(kernel)
end

