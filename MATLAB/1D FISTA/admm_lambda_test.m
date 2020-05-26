clear all
close all

% Generate example data
N = 200;
numSpots = 10;
b = zeros(N,1);
for i = 1:numSpots
   b = b + 10*rand*gaussian_basis_1D(N, N*rand, 30*rand+5);
end
b = b + 10*randn(N,1)/50;

figure(1)
plot(b)
xlabel('\theta')
ylabel('intensity')

%% Fit data

% Data Parameters
P.num_theta = size(b,1);

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 20;
P.var_theta = linspace(1/2,50,P.num_var_t).^2;

% Zero padding and mask
zPad = [0,0];
zMask = [];

% admm params
params.lambda = 0.01;
params.rho = 1;
params.isNonnegative = 1;

params.stoppingCriterion = 1;
params.maxIter = 800;
params.tolerance = 1e-8;

params.zeroPad = zPad;
params.zeroMask = zMask;

params.plotProgress = 0;

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P,zPad);
end

% Initial x
x_init = zeros(size(A0ft_stack));
for i = 1:P.num_var_t
    x_init(:,i) = b/P.num_var_t;
end

lambda_vals = logspace(-3,1,50);
M = numel(lambda_vals);
l0_norm = zeros(M,1);

for i = 1:M
    params.lambda = lambda_vals(i);
    [x_hat,err,obj] = convADMM_LASSO_1D(A0ft_stack,b,x_init,params);
    l0_norm(i) = sum(x_hat(:) > 0);
end

figure(5)
semilogx(lambda_vals,l0_norm,'o-')
xlabel('\lambda')
ylabel('l_0 norm')

b_hat = Ax_ft_1D(A0ft_stack,x_hat);

% plot 
figure(1)
hold on
plot(b_hat)

figure(2)
hold on
plot(obj)
xlabel('iteration')
ylabel('objective')

figure(3)
vdf = sum(x_hat,1)/sum(x_hat(:));
plot(vdf)