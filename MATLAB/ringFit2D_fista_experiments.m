%% Experiment setup
% Data I/O directories
data_dir = '/data/dbanco02/matlab_images/';
results_dir = '/data/dbanco02/fista_fit_results/';

% Ring sampling parameters
ring_width = 30;
num_theta= 2048;
num_rad = 2*ring_width;
dtheta = 2*pi/num_theta;
drad = 1;

% Basis function variance parameters
num_var_t = 15;
num_var_r = 10;
var_theta = linspace(dtheta,pi/32,num_var_t).^2;
var_rad   = linspace(drad,  3,       num_var_r).^2;

% Generate unshifted basis function matrices
A0ft_stack = unshifted_basis_matrix_ft_stack_norm(var_theta,var_rad,dtheta,drad,num_theta,num_rad);
A0_stack = unshifted_basis_matrix_stack_norm(var_theta,var_rad,dtheta,drad,num_theta,num_rad);

% FISTA parameters
params.stoppingCriterion = 1;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;

%% Begin experiments
for step = 1:5
    parfor img = 1:205
        disp(['step ', num2str(step),', img ' num2str(img)])
        
        % Load polar_image
        im_struc = load([data_dir,... 
        'polar_image_al7075_load_',...
        num2str(step-1), '_img_',...
        num2str(img-1), '.mat']);

        % FISTA with backtracking
        [x_hat, err, obj, l_0]  = FISTA_Circulant(A0ft_stack,im_struc.polar_image,params);   
        
        % Save result
        dir = [results_dir,...
               'fista_out_load_',...
               num2str(step-1), '_img_',...
               num2str(img-1), '.mat'];

        % Use function to save to make this valid parfor loop
        save_vars(dir,x_hat,err,obj,l_0,params);
    end
end
