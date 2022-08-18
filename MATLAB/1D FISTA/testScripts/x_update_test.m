%% Parameter selection
clear all
close all
P.set = 1;

dataset = 'D:\CHESS_data\simulated_two_spot_1D_gnoise4_nonorm_3';
output_dir = 'D:\CHESS_data\simulated_two_spot_1D_gnosie4_nonorm3_indep11';

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

% params
params.lambda = 0.01;
params.rho = 1;
params.tau = 1.1;
params.mu = 2;
params.adaptRho = 1;
params.alpha = 1.8;

params.beta = 1.01;
params.L = 1;
params.isNonnegative = 1;
params.noBacktrack = 0;

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 1600;
params.maxGradIters = 800;
params.gradTolerance = 1e-12;
params.tolerance = 1e-12;

params.numIms = 20;
params.imageNum = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;

params.plotProgress = 0;
params.verbose = 1;
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
image_num = 3;
im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));

% Zero pad image
b = im_data.polar_vector;

% Initial solution
x_init = zeros(size(A0ft_stack));
for i = 1:P.num_var_t
    x_init(:,i) = b/P.num_var_t;
end
% yk = x_init;
% vk = zeros(size(x_init));
% bnorm = norm(b(:));
% x_hat1 = circulantLinSolve( A0ft_stack/bnorm,b/bnorm,yk,vk,params );
% [x_hat2,gradIters,params,...
%      errX,l1_normX] = convGradDescent( A0ft_stack,x_init,b,yk,vk,params );
% x_hat1(x_hat1<0) = 0;
% bb1 = Ax_ft_1D(A0ft_stack,x_hat1);
% bb2 = Ax_ft_1D(A0ft_stack,x_hat2);
% 
% % plot 
% figure(1)
% hold on
% plot(bb1)
% plot(bb2)
% plot(b)
% legend('ISM','GD','Data')
 
%% Run fits
P.params.lambda = 0.01;
[x_hat1,err,obj1] = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_init,P.params);
P.params.lambda = 0.01;
[x_hat2,err,obj2] = FISTA_Circulant_1D(A0ft_stack,b,x_init,P.params);

%% Plot results
% Threshold final solution
% ft = sum(x_hat1(:))*1e-6;
% x_hat1(x_hat1(:)<ft) = 0;
% x_hat2(x_hat2(:)<ft) = 0;

ism_l0=sum(x_hat1(:)>1e-6);
fista_l0=sum(x_hat2(:)>1e-6);

bb1 = Ax_ft_1D(A0ft_stack,x_hat1);
bb2 = Ax_ft_1D(A0ft_stack,x_hat2);

% plot 
figure(1)
hold on
plot(b,'Linewidth',1)
plot(bb1)
plot(bb2)
legend('Data',sprintf('ISM %i',ism_l0),sprintf('FISTA %i',fista_l0))
title('Nonnegative (ISM ends with .1%sum(x) threshold)')

figure(2)
hold on
plot(obj1)
plot(obj2)
ylim([0 0.1])
legend(sprintf('ISM %i',ism_l0),sprintf('FISTA %i',fista_l0))
title('Objective')

% figure(3)
% subplot(2,1,1)
% imagesc(x_hat1)
% colorbar()
% title('x_{ism}')
% subplot(2,1,2)
% imagesc(x_hat2)
% colorbar()
% title('x_{fista}')

