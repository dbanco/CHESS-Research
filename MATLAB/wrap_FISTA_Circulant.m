function wrap_FISTA_Circulant(datadir,load_step,img_num,outputdir)

%wrap_FISTA_Circulant 

%% A0ft_stack
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

%% params
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;

%% load polar image

load_step
img_num
str1 = sprintf('%i',load_step);
str2 = sprintf('%i',img_num);
fileName = ['polar_image_',str1,'_',str2,'.mat'];
fileDir = fullfile(datadir,fileName);
load(fileDir)

%% call function
[x_hat, err, obj, l_0] = FISTA_Circulant(A0ft_stack,polar_image,params);

%% save output
save(fullfile(outputdir,sprintf('fista_fit_%i_%i.mat',load_step,img_num)),'x_hat','err')
end

