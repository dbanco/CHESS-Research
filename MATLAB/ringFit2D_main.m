%% Problem parameters
% Data I/O directories
data_dir = 'D:\CHESS_data\al7075_311_polar\';
results_dir = 'D:\CHESS_results\test\';

% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  2,       P.num_var_r).^2;

% Generate unshifted basis function matrices
P.betap = P.dtheta*P.drad;
P.weight = 1;
P.alphap = 10;

A0ft_stack = unshifted_basis_matrix_ft_stack(P);
A0_stack = unshifted_basis_matrix_stack(P);

% FISTA parameters
params.stoppingCriterion = 3;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 1000;
params.isNonnegative = 1;
P.params = params;

P.load_step = 4;
P.img = 35;

% Load polar_image
load([data_dir,... 
'polar_image_',...
num2str(P.load_step),'_',...
num2str(P.img), '.mat']);

% Reduce image 
test_im = polar_image(:,1:100);
test_A0ft_stack = A0ft_stack(:,1:100,:,:);

x_init = ones(size(test_A0ft_stack));
%% FISTA with backtracking
[x_hat, err, obj, l_0]  = FISTA_Circulant(test_A0ft_stack,test_im,x_init,params);

err(end)
%% Plot stuff
figure
subplot(3,1,1)
semilogy(err)
subplot(3,1,2)
semilogy(obj)
subplot(3,1,3)
semilogy(l_0)



%  save(sprintf('fista_fit_%i_%i.mat','x_hat','err','polar_image','P',...
%       P.load_step,P.img))