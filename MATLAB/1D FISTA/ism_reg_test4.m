%% Parameter selection
clear all
close all
P.set = 1;

dataset = 'D:\CHESS_data\simulated_two_spot_1D_gnoise4_nonorm_3';

num_ims = 10;
dset_name = 'gnoise4_nonorm';

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta = size(polar_vector,1);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 20;
P.var_theta = linspace(P.dtheta/2,50,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% lasso/admm params
params.lambda = 0.01;
params.rho = 1;
params.tau = 1.05;
params.mu = 2;
params.adaptRho = 1;
params.alpha = 1.8;

params.maxIter = 300;
params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.tolerance = 1e-10;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.plotProgress = 0;
params.verbose = 0;

params.maxIterReg = 1000;
params.lambda2 = 0.01;
params.rho2 = 4;
params.verbose = 0;

P.params = params;
   
baseFileName = 'fista_fit_%i_%i.mat';

% Lambda values
lambda_vals = logspace(-4,1,30);
N = numel(lambda_vals);
P.lambda_values = lambda_vals;

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
A0 = unshifted_basis_vector_stack_norm2_zpad(P);

xs = cell(num_ims,1);
bns = cell(num_ims,1);
bbs = cell(num_ims,1);

%% Run grid search ISM
x_hats = cell(3,1);
for image_num = 1:3
    im_data = load(fullfile([dataset],[prefix,'_',num2str(6*image_num),'.mat']));
    % Zero pad image
    b = im_data.polar_vector;

    % Initial solution
    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        x_init(:,i) = b/P.num_var_t;
    end

    [x_hat,~,~] = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_init,params);
    x_hats{image_num} = x_hat;
    
%     bb = Ax_ft_1D(A0ft_stack,x_hat);
%     figure(1)
%     subplot(3,1,image_num)
%     plot(b)
%     hold on
%     plot(bb)
%     title(sprintf('Diffraction image: t = %i',image_num))
    
end

%% Regularized vs unregularized comparison
% im_data = load(fullfile(dataset,[prefix,'_',num2str(2),'.mat']));
% Zero pad image
% b = im_data.polar_vector;

lam2_vals = logspace(-3,0,20);

x_n = {x_hats{1}, x_hats{3}};
params.maxIter = 700;
params.lambda = 0.1;
[x_hat1,~,~] = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_hats{2},params);
bb1 = Ax_ft_1D(A0ft_stack,x_hat1);
tvx_indep = 0;
for i = [1,3]
    tvx_indep = tvx_indep + sum(abs(x_hat1(:) - x_hats{i}(:)));
end
params.maxIterReg = 700;
params.lambda = 0.1;
params.lambda2 = 1;
params.rho2 = 4;
params.verbose = 0;

tvx_coupled = zeros(numel(lam2_vals),1);
for j = 1:numel(lam2_vals)
    params.lambda2 = lam2_vals(j);
    [x_hat2, err, obj] = convADMM_LASSO_Sherman_TVx_1D2(A0ft_stack,b,x_hats{2},x_n,params);
%     figure(2)
%     subplot(2,5,j)
%     bb2 = Ax_ft_1D(A0ft_stack,x_hat2);
%     hold on
%     plot(b)
%     plot(bb2)
    for i = [1,3]
        tvx_coupled(j) = tvx_coupled(j) + sum(abs(x_hat2(:) - x_hats{i}(:)));
    end
end

figure(1)
semilogx(lam2_vals(1),tvx_indep,'o')
hold on
semilogx(lam2_vals,tvx_coupled,'o-')
legend(sprintf( 'indep') ,...
       sprintf( 'coupled' ))
xlabel('\lambda_{TV}')
ylabel('TV')
   
%% Single test
params.plotProgress = 0;
params.verbose = 1;
params.lambda = 0;
params.lambda2 = 0;
[x_hat2, err, obj, l1_norm,tv_pen] = convADMM_LASSO_Sherman_TVx_1D2(A0ft_stack,b,x_hats{2},x_n,params);
bb2 = Ax_ft_1D(A0ft_stack,x_hat2);
figure(4)
hold on
plot(bb2)
plot(b)

figure(5)
subplot(2,2,1)
plot(obj)
title('Objective')

subplot(2,2,2)
plot(err)
title('Error')

subplot(2,2,3)
plot(l1_norm)
title('l_1 norm')

subplot(2,2,4)
plot(tv_pen)
title('TV')