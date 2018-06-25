%% Fixed Parameters
% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 100;
P.num_rad = 21;
P.num_omega = 40;

P.dtheta = 2*pi/2048;
P.drad = 1;
P.domega = 2*pi/P.num_omega;

% Basis function variance parameters
P.num_var_t = 10;
P.num_var_r = 10;
P.num_var_w = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  2,    P.num_var_r).^2;
P.var_omega = linspace(P.domega,pi/64,P.num_var_w).^2;

% fista params
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 100;
params.beta = 1.2;
params.maxIter = 20;
params.isNonnegative = 1;
P.params = params;

load('synth_3d_image.mat')

%% Generate unshifted basis function matrices
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2_3D(P);

% Initialize solution
[m,n,l,t,r,w] = size(A0ft_stack);
x_init = ones(m,n,l,t,r,w);

% Reduce image 
test_im = B(10:30,1:100,1:40);
%test_A0ft_stack = A0ft_stack(10:30,1:100,1:40);

%% FISTA with backtracking
[x_hat, err, obj, l_0]  = FISTA_Circulant_3D(A0ft_stack,test_im,x_init,params);   

save(sprintf('fista_fit_%i_%i.mat','x_hat','err','polar_image','P',...
      P.load_step,P.img))