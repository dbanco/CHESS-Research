%% Parameter selection
clear all
close all
P.set = 1;

dataset = 'E:\CHESS_data\simulated_two_spot_1D_gnoise4_nonorm_3';
output_dir = 'E:\CHESS_data\simulated_two_spot_1D_gnosie4_nonorm3_indep11';

num_ims = 10;
dset_name = 'gnoise4_nonorm';

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
load(fullfile([dataset],[prefix,'_1.mat']));
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

% ADMM params
params.lambda = 0.03;
params.rho = 0.5;
params.beta = 1.01;
params.maxIter = 800;
params.maxGradIters = 60;
params.gradTolerance = 1e-2;
params.isNonnegative = 1;
params.noBacktrack = 0;

params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 300;

0;
params.t_k = 1;

params.numIms = 20;
params.imageNum = 2;
params.zeroPad = zPad;
params.zeroMask = zMask;

params.plotProgress = 0;

P.params = params;
   
baseFileName = 'fista_fit_%i_%i.mat';

% Lambda values
lambda_vals = logspace(-3,1,30);
N = numel(lambda_vals);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
A0 = unshifted_basis_vector_stack_norm2_zpad(P);

xs = cell(num_ims,1);
bns = cell(num_ims,1);
bbs = cell(num_ims,1);

%% Run grid search GD
x_hats = cell(3,1);
for image_num = 1:3
    im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
    % Zero pad image
    b = im_data.polar_vector;


    % Initial solution
    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        x_init(:,i) = b/P.num_var_t;
    end

    [x_hat,err,obj] = convADMM_LASSO_GD_1D(A0ft_stack,b,x_init,P.params);
    x_hats{image_num} = x_hat;
    
    figure(1)
    subplot(3,1,image_num)
    plot(b)
    hold on
    plot(Ax_ft_1D(A0ft_stack,x_hat))
    title(num2str(sum(x_hat(:)>0)))
    
    figure(11)
    subplot(3,1,image_num)
    plot(obj)
end

%% Regularized
% ADMM params
params.lambda1 = 0.03;
params.lambda2 = 0.03;
params.rho1 = 0.5;
params.rho2 = 0.5;
params.beta = 1.01;
params.maxIter = 800;
params.maxGradIters = 60;
params.isNonnegative = 1;
params.noBacktrack = 0;

params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 300;
params.t_k = 1;

params.numIms = 20;
params.imageNum = 2;
params.zeroPad = zPad;
params.zeroMask = zMask;

% Load data
im_data = load(fullfile([dataset],[prefix,'_',num2str(2),'.mat']));
% Zero pad image
b = im_data.polar_vector;
x_neighbors = {x_hats{1}, x_hats{3}};


[x_hat, err, obj] = convADMM_LASSO_TVx_1D(A0ft_stack,b,x_hats{2},x_neighbors,params);

bb = Ax_ft_1D(A0ft_stack,x_hat);

figure(1)
subplot(3,1,2)
plot(bb)
sum(x_hat(:)>0)

figure(12)
plot(obj)
