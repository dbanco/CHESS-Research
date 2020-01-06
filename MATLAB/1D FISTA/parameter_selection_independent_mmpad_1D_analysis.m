%% Parameter selection
clear all
close all
P.set = 1;
% dataset = '/cluster/home/dbanco02/mmpad_polar/ring1_zero/';
% output_dir = '/cluster/shared/dbanco02/mmpad_1D_indep_param_1/';

dataset = 'D:\MMPAD_data\ring1_zero\';
output_dir = 'D:\CHESS_data\mmpad_1D_indep_param_3\';

mkdir(output_dir)
num_ims = 500;

% Universal Parameters
% Ring sampling parameters
prefix = 'mmpad_img';
load([dataset,prefix,'_1.mat']);
P.num_theta = size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r).^2;
% Zero padding and mask\
zPad = [0,0];
zMask = [129:133];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1;
params.t_k = 1;
params.lambda = 0.08;
params.wLam = 25;
params.beta = 1.2;
params.maxIter = 800;
params.maxIterReg = 800;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;
   
baseFileName = 'fista_fit_%i_%i.mat';

% Lambda values
lambda_vals = logspace(-3,1,30);

N = numel(lambda_vals);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
end
A0 = unshifted_basis_vector_stack_norm2(P);

% Select noise level
% obj_gcv = zeros(N,num_ims);
err_select = zeros(N,num_ims);

for i = 1:N
    fprintf('%i of %i\n',i,N)
    for j = 1:num_ims
        load(fullfile(output_dir,sprintf(baseFileName,i,j)))
        err_select(i,j) = err(end);
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
noise_eta = noise_est - 0.15;

% some other method
% noise_eta = max(noise_est)*norm_ratio;

discrep_crit = abs(err_select'-repmat(noise_eta,1,N));
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
    plot(fit)
    
    im_ind = im_ind + 1;
end

figure(1)
hold on
plot(awmv_az,'o')


%% Plot resulting AWMV

figure(1)
hold on
plot(awmv_az)
legend('Fit Discrep','Truth','Fit Single \lambda','Location','best')
ylabel('AWMV')
xlabel('t')

%% Plot basis functions
figure(5)
for i = 1:15
   kernel = shift1D(A0(:,i),round(cN/2));
   hold on
   plot(kernel)
end

%% Plot parameter search curves
load([param_dir,'param_search_discrep.mat'])
% lambda_vals = logspace(-3,0,500)
close all

avg_err = squeeze(mean(err_gcv,1));
avg_obj = squeeze(mean(obj_gcv,1));
std_err = squeeze(std(err_gcv,1));

figure(1)
semilogx(lambda_vals,err_gcv,'o-')
title('Average error')
ylabel('Relative Error')
xlabel('Lambda')

figure(2)
semilogx(lambda_vals,obj_gcv,'o-')
title('Average objective')
ylabel('Relative Error')
xlabel('Lambda')

figure(3)
semilogx(lambda_vals,std_err,'o-')
title('Std error')
ylabel('Relative Error')
xlabel('Lambda')

%%
figure(3)
plot(res_term,l1_term)
xlabel('Residual')
ylabel('L-1')
title('L-curve')

figure(4)
for i = 1:num_ims
    subplot(5,2,i)
    hold on
    plot(log(res_term(:,i)),log(l1_term(:,i)),'-o')
    xlabel('log-Residual')
    ylabel('log-L-1')
end

figure(5)
loglog(res_term(:,image_num),l1_term(:,image_num),'o-')
xlabel('Residual')
ylabel('L-1')

slopes = zeros(N-1,num_ims);
for i = 1:num_ims
    xx = log(res_term(:,i));
    yy = log(l1_term(:,i));
    slope = diff(yy)./diff(xx);
    slopes(:,i) = slope;
end
    
 figure(6)
for i = 1:num_ims
    subplot(5,2,i)
    hold on
    curvature = diff(slopes(2:end,i));
    yi = max(curvature);
    xi = find(curvature==yi);
    plot(curvature)
    hold on
    plot(xi,yi,'-o')
    xlabel('index')
    ylabel('slopes')
    
end

% figure(6)
% plot(res_term(:,image_num),l1_term(:,image_num),'o-')
% xlabel('Residual')
% ylabel('L-1')


% X = [res_term(:,noise_val),l1_term(:,noise_val)];
% [L,R,k] = curvature(X);
% figure(6)
% plot(k2)