%% Parameter selection
clear all
close all
P.set = 1;
dataset = ['/cluster/home/dbanco02/mmpad_polar/ring1_zero/'];
output_dir = '/cluster/shared/dbanco02/mmpad_1D_indep_param_4/';
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
zMask = [];

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

%% Run grid search
for image_num = 1:num_ims
    im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
    % Zero pad image
    b = squeeze(sum(im_data.polar_image,1));
    % Scale image by 2-norm
    bn = b/norm(b(:));
    
    % Initial solution
    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        x_init(:,i) = bn/(P.num_var_t);
    end 
 
    parfor ii = 1:N
        P_new = P;
        P_new.params.lambda = lambda_vals(ii);
        P_new.set = ii;
        [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D(A0ft_stack,bn',x_init,P_new.params);
        save_output(output_dir,baseFileName,x_hat,err,obj,im_data.polar_image,P_new,image_num)
    end
end


function save_output(output_dir,baseFileName,x_hat,err,obj,polar_image,P,image_num)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','obj','polar_image','P')
end


% %% Select lambda values
% load([param_dir,'param_search_mmpad_subset.mat']) 
% lambda_indices = zeros(num_ims,1);
% noise_est = zeros(num_ims,1);
% norm_data = zeros(num_ims,1);
% noise_eta = zeros(num_ims,1);
% for image_num = 1:num_ims
%     im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
%     % Zero pad image
%     b = im_data.polar_image;
%     % Scale image by 2-norm
%     bn = b/norm(b(:));
%     
%     [rN,cN] = size(b);
%     kernel = shift2D(A0(:,:,1,1),round(rN/2),round(cN/2));
%     kernel = kernel./sum(kernel(:));
%     kernel = kernel(1:35,:);
%     bn_hat = conv2(bn,kernel,'same');
%     noise_est(image_num) = norm(bn-bn_hat);
%     norm_data(image_num) = norm(b);
% end
% % norm scaling method
% noise_eta = norm_data./max(norm_data).*max(noise_est);
% % purely residual based method
% % noise_eta = noise_est;
% 
% discrep_crit = abs(err_gcv'-repmat(noise_eta,1,N));
% [lambda_indices,~] = find(discrep_crit' == min(discrep_crit'));
% param_select = lambda_vals(lambda_indices);
% figure(2)
% plot(param_select)
% save([param_dir,'lambda_select_mmpad_subset1.mat'],'param_select')
% %% Compute true awmv
% load(['D:\CHESS_data\simulated_two_spot_1D_noise2_12\synth_data.mat'])
% truth_awmv_az = zeros(num_ims,1);
% 
% for i = 1:num_ims
%     sample = synth_sample{i};
%     az_awmv = sample.std_theta'*sample.amplitudes/sum(sample.amplitudes(:));
%     truth_awmv_az(i) = az_awmv;
% end
% 
% %% Fit with selected parameters and plot
% output_dir = 'D:\CHESS_data\noise_1D_indep_param_8';
% load([param_dir,'lambda_select_noise8.mat'])
% figure(222)
% [ha1, pos1] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
% awmv_az = zeros(num_ims,1);
% for image_num = 1:num_ims
%     im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
%     % Zero pad image
%     b = im_data.polar_vector;
%     % Scale image by 2-norm
%     bn = b/norm(b(:));
%     
%     % Initial solution
%     x_init = zeros(size(A0ft_stack));
%     for i = 1:P.num_var_t
%         x_init(:,i) = bn/P.num_var_t;
%     end
%     
%     P.params.lambda = param_select(image_num);
%     [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D_Poisson(A0ft_stack,bn,x_init,P.params);
%     fit = Ax_ft_1D(A0ft_stack,x_hat);
%     var_signal = squeeze(sum(x_hat,1));
%     var_sum = sum(var_signal(:));
%     awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*var_signal(:))/var_sum;
%     polar_image = polar_vector;
%     
%     save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','polar_image','P')
%     
%     % Plot
%     axes(ha1(image_num))
%     hold on
%     
%     plot(bn)
%     plot(fit)
% end
% 
% figure(1)
% hold on
% plot(awmv_az)
% plot(truth_awmv_az)
% 
% %% Fit single parameter
% figure(333)
% [ha2, pos2] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
% awmv_az = zeros(num_ims,1);
% P.params.lambda = mean(param_select);
% for image_num = 1:num_ims
%     im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
%     % Zero pad image
%     b = im_data.polar_vector;
%     % Scale image by 2-norm
%     bn = b/norm(b(:));
%     
%     % Initial solution
%     x_init = zeros(size(A0ft_stack));
%     for i = 1:P.num_var_t
%         x_init(:,i) = bn/P.num_var_t;
%     end
%     
%     [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D_Quadratic(A0ft_stack,bn,x_init,P.params);
%     fit = Ax_ft_1D(A0ft_stack,x_hat);
%     var_signal = squeeze(sum(x_hat,1));
%     var_sum = sum(var_signal(:));
%     awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*var_signal(:))/var_sum;
%     
%     % Plot
%     axes(ha2(image_num))
%     hold on
%     
%     plot(bn)
%     plot(fit)
%     
% end
% 
% figure(1)
% hold on
% plot(awmv_az)
% legend('Fit Discrep','Truth','Fit Single \lambda','Location','best')
% ylabel('AWMV')
% xlabel('t')
% 
% %% Plot parameter search curves
% load([param_dir,'param_search_discrep.mat'])
% % lambda_vals = logspace(-3,0,500)
% close all
% 
% avg_err = squeeze(mean(err_gcv,1));
% avg_obj = squeeze(mean(obj_gcv,1));
% std_err = squeeze(std(err_gcv,1));
% 
% figure(1)
% semilogx(lambda_vals,err_gcv,'o-')
% title('Average error')
% ylabel('Relative Error')
% xlabel('Lambda')
% 
% figure(2)
% semilogx(lambda_vals,obj_gcv,'o-')
% title('Average objective')
% ylabel('Relative Error')
% xlabel('Lambda')
% 
% figure(3)
% semilogx(lambda_vals,std_err,'o-')
% title('Std error')
% ylabel('Relative Error')
% xlabel('Lambda')
% 
% %%
% figure(3)
% plot(res_term,l1_term)
% xlabel('Residual')
% ylabel('L-1')
% title('L-curve')
% 
% figure(4)
% for i = 1:num_ims
%     subplot(5,2,i)
%     hold on
%     plot(log(res_term(:,i)),log(l1_term(:,i)),'-o')
%     xlabel('log-Residual')
%     ylabel('log-L-1')
% end
% 
% figure(5)
% loglog(res_term(:,image_num),l1_term(:,image_num),'o-')
% xlabel('Residual')
% ylabel('L-1')
% 
% slopes = zeros(N-1,num_ims);
% for i = 1:num_ims
%     xx = log(res_term(:,i));
%     yy = log(l1_term(:,i));
%     slope = diff(yy)./diff(xx);
%     slopes(:,i) = slope;
% end
%     
%  figure(6)
% for i = 1:num_ims
%     subplot(5,2,i)
%     hold on
%     curvature = diff(slopes(2:end,i));
%     yi = max(curvature);
%     xi = find(curvature==yi);
%     plot(curvature)
%     hold on
%     plot(xi,yi,'-o')
%     xlabel('index')
%     ylabel('slopes')
%     
% end
% 
% % figure(6)
% % plot(res_term(:,image_num),l1_term(:,image_num),'o-')
% % xlabel('Residual')
% % ylabel('L-1')
% 
% 
% % X = [res_term(:,noise_val),l1_term(:,noise_val)];
% % [L,R,k] = curvature(X);
% % figure(6)
% % plot(k2)