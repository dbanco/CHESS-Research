%% Parameter selection
clear all
close all
P.set = 1;
datadir = 'D:\CHESS_data\';
dataset = 'D:\CHESS_data\simulated_two_spot_1D_noise2';
num_ims = 10;

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
load(fullfile([dataset,'_1'],[prefix,'_1.mat']));
P.num_theta = size(polar_vector,1);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 15;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1;
params.t_k = 1;
params.lambda = 0;
params.wLam = 25;
params.gamma = 1;
params.maxIterReg = 800;
params.beta = 1.2;
params.maxIter = 800;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;

P.params = params;

% coupled params
load('lambda_select1.mat')
Pc.lambda_values = param_select;
Pc.wLam = 25;
Pc.gamma = 1;
Pc.maxIterReg = 800;
Pc.num_outer_iters = 10;
Pc.baseFileName = 'fista_fit_%i_%i.mat';
Pc.num_ims = num_ims;
Pc.prefix = 'polar_vector';
Pc.dataset = [dataset,'_2'];
% Gamma values
gamma_vals = [0.05,0.1, 0.2]; 
N = numel(gamma_vals);

% Select noise level
dir_num = sprintf('_%i',1);
obj_gcv = zeros(N,num_ims);
err_gcv = zeros(N,num_ims);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
end
A0 = unshifted_basis_vector_stack_norm2(P);
%% Run grid search
for i = 1:N
    Pc.init_dir = [datadir,'noise_1D_coupled_init',num2str(i)];
    Pc.output_dirA = [datadir,'noise_1D_coupled_',num2str(i),'a'];
    Pc.output_dirB = [datadir,'noise_1D_coupled_',num2str(i),'b'];
    mkdir(Pc.init_dir)
    mkdir(Pc.output_dirA)
    mkdir(Pc.output_dirB)
    Pc.gamma = gamma_vals(i);
    runCoupledFISTA_1D(P,Pc)
end
% save('param_search_coupled_discrep_1.mat','err_gcv','obj_gcv','P','lambda_vals')

%% Select lambda values
err_discrep = zeros(num_ims,N);
noise_eta = zeros(num_ims,N);
for i = 1:N
    for image_num = 1:num_ims
        im_data = load([datadir,'noise_1D_coupled_',num2str(i),'a\',sprintf(Pc.baseFileName,1,image_num)]);
        err_discrep(image_num,i) = im_data.err(end-1);
        b = im_data.polar_image;
        % Scale image by 2-norm
        bn = b/norm(b(:));

    %     kernel = [0.0003 0.1329 0.9822 0.1329 0.0003];
        kernel = [0.0001 0.0003 0.0012 0.0042 0.0127 0.0329,...
                  0.0740 0.1434 0.2399 0.3466 0.4322 0.4652];
        kernel = [kernel,fliplr(kernel(1:end-1))];
        kernel = kernel./sum(kernel(:));

        bn_hat = conv(bn,kernel,'same');
        noise_eta(image_num,i) = norm(bn-bn_hat);
    end
end
discrep_crit = abs(squeeze(sum(err_discrep-noise_eta,1)))
gamma_index = find(discrep_crit == min(discrep_crit));


%% Compute true awmv
load(['D:\CHESS_data\simulated_two_spot_1D_noise2_12\synth_data.mat'])
truth_awmv_az = zeros(num_ims,1);

for i = 1:num_ims
    sample = synth_sample{i};
    az_awmv = sample.std_theta'*sample.amplitudes/sum(sample.amplitudes(:));
    truth_awmv_az(i) = az_awmv;
end

%% Fit with selected parameters and plot
figure(222)
[ha1, pos1] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az = zeros(num_ims,1);
for image_num = 1:num_ims
    im_data = load(fullfile([datadir,dir_num],[prefix,'_',num2str(image_num),'.mat']));
    % Zero pad image
    b = im_data.polar_vector;
    % Scale image by 2-norm
    bn = b/norm(b(:));
    
    % Initial solution
    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        x_init(:,i) = bn/P.num_var_t;
    end
    
    P.params.lambda = param_select(image_num);
    [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D_Quadratic(A0ft_stack,bn,x_init,P.params);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    var_signal = squeeze(sum(x_hat,1));
    var_sum = sum(var_signal(:));
    awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*var_signal(:))/var_sum;
    
    % Plot
    axes(ha1(image_num))
    hold on
    
    plot(bn)
    plot(fit)
end

figure(1)
hold on
plot(awmv_az)
plot(truth_awmv_az)

%% Fit single parameter
figure(333)
[ha2, pos2] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az = zeros(num_ims,1);
P.params.lambda = mean(param_select);
for image_num = 1:num_ims
    im_data = load(fullfile([dataset,dir_num],[prefix,'_',num2str(image_num),'.mat']));
    % Zero pad image
    b = im_data.polar_vector;
    % Scale image by 2-norm
    bn = b/norm(b(:));
    
    % Initial solution
    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        x_init(:,i) = bn/P.num_var_t;
    end
    
    [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D_Quadratic(A0ft_stack,bn,x_init,P.params);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    var_signal = squeeze(sum(x_hat,1));
    var_sum = sum(var_signal(:));
    awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*var_signal(:))/var_sum;
    
    % Plot
    axes(ha2(image_num))
    hold on
    
    plot(bn)
    plot(fit)
    
end

figure(1)
hold on
plot(awmv_az)
legend('Fit Discrep','Truth','Fit Single \lambda','Location','best')
ylabel('AWMV')
xlabel('t')

%% Plot parameter search curves
load('param_search_discrep.mat')
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