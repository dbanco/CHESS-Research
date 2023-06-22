%% Generate example data
% Gaussian peak functions with Poisson statistics
close all

T = 1;
N = 100;
[b,~,awmv_true] = generateExampleData2D(N,N,T);

figure(1)
subplot(3,1,1)
imagesc(b(:,:,1))
title('data')

%% Define parameters
P.N1 = size(b,1);
P.N2 = size(b,2);

% Basis function variance parameters
P.K1 = 10;
P.K2 = 15;
P.K = P.K1*P.K2;
P.sigma1 = linspace(0.5,  10,   P.K1);
P.sigma2 = linspace(0.5,  10,   P.K2);
P.basis = 'norm2';
P.mu1 = round(P.N1/2);
P.mu2 = round(P.N2/2);

% fista params
params.lambda = 0.05; % sparsity penalty
params.L = 1000;  %
params.beta = 2; %
params.noBacktrack = 0;

params.isNonnegative = 1; % flag to enforce nonnegativity

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 10;
params.tolerance = 1e-3;

params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
params.verbose = 1;      % flag to print objective values at each iteration 
P.params = params;

% gcmp params
opt = struct();
opt.delta = 1e-3;
opt.epsilon = 1;
opt.maxLayers = 1;
opt.NonnegativeDict = 1;
opt.NonNegCoef = 1;
opt.DictFilterSizes = ones(2,P.K1*P.K2);

j = 1;
for k1 = 1:P.K1
    for k2 = 1:P.K2
        opt.DictFilterSizes(1,j) = ceil(5*(P.sigma1(k1)));
        opt.DictFilterSizes(2,j) = ceil(5*(P.sigma2(k2)));
        j = j + 1;
    end
end

% Construct dictionary
D = dictionary2D(P);
Df = fft2(D);
% Initialize solution
x_init = zeros(P.N1,P.N2,P.K);

%% Solve FISTA
X1 = FISTA_Circulant_cpu(Df,b,x_init,params);
% Comptue reconstruction of image
Xf1 = fft2(X1);
b_hat1 = ifft2(Ax_cpu(Df,Xf1),'symmetric');


%% Solve GCMP
X2 = solveGCMP(Df,b,opt);
% Comptue reconstruction of image
Xf2 = fft2(X2);
b_hat2 = ifft2(Ax_cpu(Df,Xf2),'symmetric');


%% Display results
fprintf('FISTA error: %0.2f \n',norm(b(:)-b_hat1(:))/norm(b(:)) )
fprintf('FISTA sparsity: %i nonzeros \n',sum(X1(:)>0))
fprintf('GCMP error: %0.2f \n',norm(b(:)-b_hat2(:))/norm(b(:)) )
fprintf('GCMP sparsity: %i nonzeros \n',sum(X2(:)>0))

subplot(3,1,2)
imagesc(b_hat1)
title('FISTA solution')

subplot(3,1,3)
imagesc(b_hat2)
title('GCMP Solution')






