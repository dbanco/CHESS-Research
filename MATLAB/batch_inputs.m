% Data directory
datadir = fullfile('cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar');

% Job directory
jobdir = fullfile(datadir,'job_al7075_311');
mkdir(jobdir)

% Output directory
outputdir = fullfile(datadir,'al7075_311_polar_fit');
mkdir(outputdir)

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
A0ft_stack = unshifted_basis_matrix_ft_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad);
A0_stack = unshifted_basis_matrix_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad);
save(fullfile(jobdir,'basis.mat'),'A0ft_stack')

% FISTA parameters
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;
save(fullfile(jobdir,'params.mat'),'params')

%% Parameters to vary
steps = 1:5;
imgs = 1:205;
k = 0;
for step = steps
    for img = img
        k = k + 1;
        save(fullfile(jobdir,['varin_',num2str(k),'.mat']),'step','img','outputdir')
    end
end
        
num_jobs = k;