%% Problem parameters
% Data I/O directories
data_dir = 'D:\CHESS_data\al7075_311_polar\';
results_dir = 'D:\CHESS_results\fista_fit_results\';

% Ring sampling parameters
ring_width = 20;
num_theta= 2048;
num_rad = 2*ring_width+1;
dtheta = 2*pi/num_theta;
drad = 1;

% Basis function variance parameters
num_var_t = 15;
num_var_r = 10;
var_theta = linspace(dtheta,pi/64,num_var_t).^2;
var_rad   = linspace(drad,  2,       num_var_r).^2;

% Generate unshifted basis function matrices
betap = 1;
weight = 1;
A0ft_stack = unshifted_basis_matrix_ft_stack_weight(var_theta,var_rad,dtheta,drad,num_theta,num_rad,betap);
A0_stack = unshifted_basis_matrix_stack_weight(var_theta,var_rad,dtheta,drad,num_theta,num_rad,betap);

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