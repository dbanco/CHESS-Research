% Load data
data_dir = 'E:\CHESS_data\al7075_311_polar';
file_name = 'polar_image_%i_%i.mat';

set = 2;
img = 30;

load(fullfile( data_dir, sprintf(file_name,set,img) ))

% Parameters
P.set = set;
P.num_rad = size(polar_image,1);
P.num_theta = size(polar_image,2);

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(1,50,P.num_var_t).^2;
P.var_rad   = linspace(1,  3,    P.num_var_r).^2;
P.basis = 'norm2';

% fista params
params.lambda = 0.08; % sparsity penalty
params.L = 75000;  %
params.beta = 2; %
params.noBacktrack = 0;

params.isNonnegative = 1; % flag to enforce nonnegativity

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 100;
params.tolerance = 1e-6;

% Zero padding and mask (just ignore this)
zPad = [0,0];
zMask = [];
params.zeroPad = zPad; % number of [rows,columns]of padding to add to data
params.zeroMask = zMask; % specifies columns known to be zero

params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
params.verbose = 1;      % flag to print objective values at each iteration 
P.params = params;


%% Setup and solve
% Construct dictionary
% A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
A0ft_stack = unshifted_basis_matrix_ft_stack(P);

% Initialize solution
x_init = zeros(size(A0ft_stack));
% b = sum(polar_image/5000,1)';
% b = b(1:2047);
b = polar_image/100;

% Solve
[x_hat,err,obj] = FISTA_Circulant(A0ft_stack,b,x_init,params);

%% View images
fit = Ax_ft_2D(A0ft_stack,x_hat);
figure(1)
subplot(2,1,1)
imagesc(b)
subplot(2,1,2)
imagesc(fit)

%% View fit 
% fit = Ax_ft_1D(A0ft_stack,x_hat);
% plot(b)
% hold on
% plot(fit)
% 
% legend('data','fit')

