% Test 1D FISTA
inDir = 'D:\CHESS_data\simulated_two_spot_1D\polar_vector';
outDir = 'D:\CHESS_data\simulated_two_spot_1D_fit';
mkdir(outDir)

P.dtheta = 1;
P.sampleDims = [71,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 100;
P.var_theta = linspace(P.dtheta/2,  32,P.num_var_t).^2;
P.num_theta = size(polar_vector,2);
    
% Zero padding and mask
zPad = [];
zMask = [];

A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);

params.stoppingCriterion = 1;
params.tolerance = 1e-10;
params.L = 1;
params.t_k = 1;
params.lambda = 0.15;
params.beta = 1.01;
params.maxIter = 5000;
params.maxIterReg = 50;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

load('D:\CHESS_data\simulated_two_spot_1D_fit\fista_fit_1_1.mat')
b = polar_vector/norm(polar_vector);
0.5*norm(b-Ax_ft_1D(A0ft_stack,x_hat))^2 + params.lambda*norm(x_hat(:),1)

truth = 2*gaussian_basis_1D_norm2( P.num_theta, 50, 10^2);
0.5*norm(polar_vector-truth)^2 + params.lambda*2

figure(3)
var_theta   = linspace(P.dtheta/2,  32,P.num_var_t);
plot(var_theta,sum(x_hat,1))