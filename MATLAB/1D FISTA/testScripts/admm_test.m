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

% fista params
params.lambda = 0.03;
params.rho = 0.5;

params.beta = 1.01;
params.L = 1;
params.isNonnegative = 0;
params.noBacktrack = 0;

params.stoppingCriterion = 2;
params.maxIter = 800;
params.maxGradIters = 200;
params.gradTolerance = 1e-1;
params.tolerance = 1e-8;

params.numIms = 20;
params.imageNum = 1;
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
image_num = 1;
im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));

% Zero pad image
b = im_data.polar_vector;

% Initial solution
x_init = zeros(size(A0ft_stack));
for i = 1:P.num_var_t
    x_init(:,i) = b/P.num_var_t;
end

[x_hat,err,obj] = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_init,P.params);

bb = Ax_ft_1D(A0ft_stack,x_hat);
xs{image_num} = x_hat;
bbs{image_num} = bb;
bns{image_num} = b;

vdfs = zeros(20,10);
xi = xs{image_num};
vdfs(:,image_num) = sum(xi,1)./sum(xi(:)); 

% plot 
figure(1)
hold on
plot(bbs{image_num})
plot(bns{image_num})

%%
% params.tolerance = 1e-6;
% A = rand(4,20); b = rand(4,1); x_in = zeros(20,1);
% [x err obj l_0] = GD_test(A,b,x_in,params);
% [A*x,b]
%% Run grid search FISTA
% 
% image_num = 1;
% im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
% % Zero pad image
% b = im_data.polar_vector;
% 
% 
% % Initial solution
% x_init = zeros(size(A0ft_stack));
% for i = 1:P.num_var_t
%     x_init(:,i) = b/P.num_var_t;
% end
% P.params.maxIter = 400;
% [x_hat,err,obj,~,~] = FISTA_Circulant_1D(A0ft_stack,b,x_init,P.params);
% 
% bb = Ax_ft_1D(A0ft_stack,x_hat);
% xs{image_num} = x_hat;
% bbs{image_num} = bb;
% bns{image_num} = b;
% 
% vdfs = zeros(20,10);
% xi = xs{image_num};
% vdfs(:,image_num) = sum(xi,1)./sum(xi(:)); 
% 
% % plot 
% 
% figure(2)
% hold on
% plot(bbs{image_num})
% plot(bns{image_num})

