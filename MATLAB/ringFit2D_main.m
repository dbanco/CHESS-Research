%% Problem parameters
% Data I/O directories
data_dir = 'D:\CHESS_data\al7075_311_polar\';
results_dir = 'D:\CHESS_results\fista_fit_results\';

% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(dtheta,pi/64,num_var_t).^2;
P.var_rad   = linspace(drad,  2,       num_var_r).^2;

% Generate unshifted basis function matrices
P.betap = P.dtheta*P.drad;
P.weight = 1;
P.alphap = 10;

A0ft_stack = unshifted_basis_matrix_ft_stack(P);
A0_stack = unshifted_basis_matrix_stack(P);

% FISTA parameters
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;

step = 4;
img = 35;
% Load polar_image
load([data_dir,... 
'polar_image_',...
num2str(step),'_',...
num2str(img), '.mat']);

%% FISTA with backtracking
[x_hat, err, obj, l_0]  = FISTA_Circulant(A0ft_stack,polar_image,params);   

save('fista_fit_wavlt_4_35.mat','x_hat','err','betap','weight','dtheta','drad','num_rad','num_theta','var_theta','var_rad','params','img','step')