%% Setup example problem
% Ring sampling parameters
ring_width = 15;
P.num_theta= 60;
P.num_rad = 2*ring_width+1;
P.dtheta = 1;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 4;
P.num_var_r = 3;
P.var_theta = linspace(P.dtheta, 8, P.num_var_t).^2;
P.var_rad   = linspace(P.drad,   3, P.num_var_r).^2;

% Create dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);

% Create spots
% Generate diffraction spot parameters
P.theta_mean1 = round(P.num_theta/4);
P.theta_mean2 = round(3*P.num_theta/4);
P.rad_mean = round(P.num_rad/2);

% Generate image
B1 = gaussian_basis_wrap_2D_norm2( P.num_theta,P.dtheta,P.theta_mean1,P.var_theta(2),...
                            P.num_rad,  P.drad,  P.rad_mean,  P.var_rad(1));
                        
B2 = gaussian_basis_wrap_2D_norm2( P.num_theta,P.dtheta,P.theta_mean2,P.var_theta(2),...
                            P.num_rad,  P.drad,  P.rad_mean,  P.var_rad(1));
B = 100*B1 + 20*B2;

%% Run GCMP
x = B;
params.epsilon = 0.1;
params.showImage = 1;
params.isNonnegative = 1;
params.zeroMask = [];

a = GCMP(A0ft_stack,x,params);