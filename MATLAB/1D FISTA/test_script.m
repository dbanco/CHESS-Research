
clear all
close all
tic
%% Generate example data
N = 201;
numSpots = 1;
b = zeros(N,1);
amplitude = 5;
mean_param = N*0.5;
var_param = 20;
for i = 1:numSpots
   b = b + amplitude*gaussian_basis_1D(N, mean_param, var_param);
end

% Add noise
b = b + 100*randn(N,1)/50;

% % Remove negative values
% b(b<0) = 0;

%% Define parameters

% Length of intensity data (theta coordinate)
P.num_theta = size(b,1); 

% Define dictionary of Gaussian basis functions
P.num_var_t = 20;   % Number of different basis functions 
P.var_theta = linspace(1/2,50,P.num_var_t).^2; % Variances of basis functions
P.basis = 'norm2';
% Zero padding and mask (just ignore this)
zPad = [0,0];
zMask = [];

% ADMM parameters
lambda1 = 1e-2;
params.lambda1 = lambda1; % sparsity penalty
params.rho1 = 0.001;  % initial ADMM


params.adaptRho = 1; % binary flag for adaptive rho
params.mu = 2;       % tolerated factor between primal/dual residual
params.tau = 1.05;   % rho update factor

params.alpha = 1.8; % over-relaxation paramter

params.isNonnegative = 1; % flag to enforce nonnegativity

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 20;
params.conjGradIter = 100;
params.tolerance = 1e-8;
params.cgEpsilon = 1e-3;

params.zeroPad = zPad; % number of [rows,columns]of padding to add to data
params.zeroMask = zMask; % specifies columns known to be zero

params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
params.verbose = 1;      % flag to print objective values at each iteration 
P.params = params;
%% Setup and solve

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
A0_stack = unshifted_basis_vector_stack_zpad(P);
% Initialize solution
x_init = zeros(size(A0ft_stack));

% Solve
[x_hat1,err,obj] = convADMM_LASSO_Sherman_1D(A0ft_stack/norm(b),b/norm(b),x_init,params);

% Compute result
b_hat = Ax_ft_1D(A0ft_stack,x_hat1);

toc
% Plot fit
figure(1)
subplot(1,2,1)
plot(b)
hold on
plot(b_hat)
xlabel('\theta')
ylabel('Intensity')
title('Data fit')
legend('data','fit')

% Plot variance distribution function
vdf = sum(x_hat1,1)/sum(x_hat1,'all');
subplot(1,2,2)
bar(vdf)
xlabel('narrow --> \sigma index --> wide')
ylabel('\Sigma x_i / \Sigma x')
title('VDF')

%% Indep CG
[x_hat2,err,obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_1D(A0ft_stack/norm(b),b/norm(b),x_init,params);
% [x_hat2,err,obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_1D(A0ft_stack, b, x_init, params);

%% Coupled CG
params.rho2 = 0;
params.lambda2 = 0;
params.lambda1 = lambda1*ones(3,1);
B = [b,b+1,b+2]; 
[X_hat3,err,obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack/norm(b),B/norm(b),zeros(N,20,3),params);
x_hat3 = X_hat3(:,:,1);

% Compute result
b_hat2 = Ax_ft_1D(A0ft_stack,x_hat2);
% Plot fit
figure(11)
subplot(1,2,1)
plot(b)
hold on
plot(b_hat2)
xlabel('\theta')
ylabel('Intensity')
title('Data fit')
legend('data','fit')

% Plot variance distribution function
vdf2 = sum(x_hat2,1)/sum(x_hat2,'all');
subplot(1,2,2)
bar(vdf2)
xlabel('narrow --> \sigma index --> wide')
ylabel('\Sigma x_i / \Sigma x')
title('VDF')

% Compute result
b_hat3 = Ax_ft_1D(A0ft_stack,x_hat3);
% Plot fit
figure(111)
subplot(1,2,1)
plot(b)
hold on
plot(b_hat3)
xlabel('\theta')
ylabel('Intensity')
title('Data fit')
legend('data','fit')

% Plot variance distribution function
vdf3 = sum(x_hat3,1)/sum(x_hat3,'all');
subplot(1,2,2)
bar(vdf3)
xlabel('narrow --> \sigma index --> wide')
ylabel('\Sigma x_i / \Sigma x')
title('VDF')
%{
Notes:

Changing the lambda1 parameter will render different results:
- lambda1 = 1 fits the peak well and not the noise (this is ideal)
- lambda1 = 0.1 fits the peak and the noise well
- lambda1 = 0.01 fits everything almost exactly
- lambda1 = 10 fits peak poorly and fits mean of noise with wide basis
            function
- lambda1 = 100 fit is worse

We can look at which basis functions were used and how that varies with
selection of lambda1:
- lambda1 = 1 
    Basis function 2 with sigma ~3.1 contributes ~40% of the signal
    Basis function 3 with sigma ~5.7 contributes ~53% of the signal
    True sigma is sqrt(20) = 4.47
    3.1*0.40 + 5.7*0.53 = 4.3
    The 2 basis functions nearest in size to the peak in the data
    contribute 93% of the signal
- lambda1 = 0.01 Biases solution towards using the smallest basis function 
    which allows the algorithm to fit everything almost exactly
- lambda1 = 100 Biases solution towards using the larger basis functions
    than are necessary because they tend to be cheaper in terms of l1-norm
    cost.

%}

