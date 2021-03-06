clear all
close all

%% Generate example data
N = 200;
numSpots = 1;
b = zeros(N,1);
amplitude = 5;
mean_param = N*0.5;
var_param = 20;
for i = 1:numSpots
   b = b + amplitude*gaussian_basis_1D(N, mean_param, var_param);
end

% Add noise
b = b + 10*randn(N,1)/50;

% Remove negative values
b(b<0) = 0;

%% Define parameters

% Length of intensity data (theta coordinate)
P.num_theta = size(b,1); 

% Define dictionary of Gaussian basis functions
P.num_var_t = 20;   % Number of different basis functions 
P.var_theta = linspace(1/2,50,P.num_var_t).^2; % Variances of basis functions

% Zero padding and mask (just ignore this)
zPad = [0,0];
zMask = [];

% ADMM parameters
params.lambda1 = 10; % sparsity penalty
params.rho1 = 0.5;  % initial ADMM


params.adaptRho = 1; % binary flag for adaptive rho
params.mu = 2;       % tolerated factor between primal/dual residual
params.tau = 1.05;   % rho update factor

params.alpha = 1.8; % over-relaxation paramter

params.isNonnegative = 1; % flag to enforce nonnegativity

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 800;
params.tolerance = 1e-6;

params.zeroPad = zPad; % number of [rows,columns]of padding to add to data
params.zeroMask = zMask; % specifies columns known to be zero

params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
params.verbose = 1;      % flag to print objective values at each iteration 
P.params = params;
%% Setup and solve

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);

% Initialize solution
x_init = zeros(size(A0ft_stack));

% Solve
[x_hat,err,obj] = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_init,params);

% Compute result
b_hat = Ax_ft_1D(A0ft_stack,x_hat);

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
vdf = sum(x_hat,1)/sum(x_hat,'all');
subplot(1,2,2)
bar(vdf)
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
