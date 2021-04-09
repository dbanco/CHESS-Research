
clear all
close all
tic
%% Generate example data
N = 121;
T = 10;
numSpots = 1;
b = zeros(N,T);
amplitude = 5;
mean_param = N*0.5;
var_param = (2:11).^2;
for t = 1:T
    for i = 1:numSpots
       b(:,t) = b(:,t) + amplitude*gaussian_basis_1D(N, mean_param, var_param(t));
    end
    % Add noise
    b(:,t) = b(:,t) + 0.1*randn(N,1);
end

data_fig = figure(1);
[ha2, ~] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
indx = 1;
for t = 1:T
    axes(ha2(indx))
    plot(b(:,t))
    ylim([0,6])
    indx = indx + 1;
end

%% Define parameters

% Length of intensity data (theta coordinate)
P.num_theta = size(b,1); 

% Define dictionary of Gaussian basis functions
P.num_var_t = 20;   % Number of different basis functions 
K = P.num_var_t;
P.var_theta = linspace(1/2,50,P.num_var_t).^2; % Variances of basis functions
P.basis = 'norm2';
% Zero padding and mask (just ignore this)
zPad = [0,0];
zMask = [];

% ADMM parameters
params.lambda1 = 5e-3; % sparsity penalty
params.rho1 = 0.001;  % initial ADMM


params.adaptRho = 1; % binary flag for adaptive rho
params.mu = 2;       % tolerated factor between primal/dual residual
params.tau = 1.05;   % rho update factor

params.alpha = 1.8; % over-relaxation paramter

params.isNonnegative = 1; % flag to enforce nonnegativity

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 100;
params.conjGradIter = 50;
params.tolerance = 1e-8;
params.cgEpsilon = 1e-5;

params.zeroPad = zPad; % number of [rows,columns]of padding to add to data
params.zeroMask = zMask; % specifies columns known to be zero

params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
params.verbose = 1;      % flag to print objective values at each iteration 
P.params = params;
%% Setup and solve

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Initialize solution
x_init = zeros(size(A0ft_stack));
b_hat = zeros(N,T);
X_hat1 = zeros(N,K,T);
for t  = 1:T
    % Solve
    bb = b(:,t)
    [x_hat1,err,obj] = convADMM_LASSO_Sherman_1D(A0ft_stack/norm(bb),bb/norm(bb),x_init,params);
    figure(8)
    plot(obj)
    % Compute result
    X_hat1(:,:,t) = x_hat1;
    b_hat(:,t) = Ax_ft_1D(A0ft_stack,x_hat1);
end

awmv1 = zeros(T,1);
var_signal = squeeze(sum(X_hat1,1));
var_sum = squeeze(sum(var_signal,1));
for t = 1:T
    awmv1(t) = sum(sqrt(P.var_theta(:)).*var_signal(:,t))/var_sum(t);
end

%% Plot over data
figure(1)
for t = 1:T
    axes(ha2(t))
    hold on
    plot(b_hat(:,t))
    ylim([0,6])
end

%% Solve using coupled algorithm
P.params.rho2 = 0.001;
P.params.lambda1 = ones(T,1)*params.lambda1;
P.params.lambda2 = 5e-2;
P.params.rho2 = 1e-2;
[X_hat2,err,obj,l1_norm,tv_penalty] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,b,X_hat1,P.params); 
figure(9)
plot(obj)
awmv2 = zeros(T,1);
var_signal = squeeze(sum(X_hat2,1));
var_sum = squeeze(sum(var_signal,1));
for t = 1:T
    awmv2(t) = sum(sqrt(P.var_theta(:)).*var_signal(:,t))/var_sum(t);
end

%% Plot over data
couple_fig = figure(2);
[ha3, ~] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
indx = 1;
for t = 1:T
    b_hat2(:,t) = Ax_ft_1D(A0ft_stack,squeeze(X_hat2(:,:,t)));
    
    axes(ha3(t))
    hold on
    plot(b(:,t))
    plot(b_hat2(:,t))
    ylim([0,6])
    indx = indx + 1;
end

figure(3)
hold on 
truth = sqrt(var_param);
plot(truth,'s-','Linewidth',2')
plot(awmv1,'o-')
plot(awmv2,'x-')
legend('truth','indep','coupled')
% toc
% % Plot fit
% figure(1)
% subplot(1,2,1)
% plot(b)
% hold on
% plot(b_hat)
% xlabel('\theta')
% ylabel('Intensity')
% title('Data fit')
% legend('data','fit')
% 
% % Plot variance distribution function
% vdf = sum(x_hat1,1)/sum(x_hat1,'all');
% subplot(1,2,2)
% bar(vdf)
% xlabel('narrow --> \sigma index --> wide')
% ylabel('\Sigma x_i / \Sigma x')
% title('VDF')

%%
% [x_hat2,err,obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_1D(A0ft_stack/norm(b),b/norm(b),x_init,params);
% % Compute result
% b_hat2 = Ax_ft_1D(A0ft_stack,x_hat2);
% % Plot fit
% figure(11)
% subplot(1,2,1)
% plot(b)
% hold on
% plot(b_hat2)
% xlabel('\theta')
% ylabel('Intensity')
% title('Data fit')
% legend('data','fit')
% 
% % Plot variance distribution function
% vdf2 = sum(x_hat2,1)/sum(x_hat2,'all');
% subplot(1,2,2)
% bar(vdf2)
% xlabel('narrow --> \sigma index --> wide')
% ylabel('\Sigma x_i / \Sigma x')
% title('VDF')
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

