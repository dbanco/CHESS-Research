center = [1025,1020];
r1 = 430;
r2 = 450;
% factor = 5;
% fname = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\ff_000277.ge2';
% onlineDir = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\277\';

fname = 'C:\Users\dpqb1\Documents\Data\c103_Feb\ff_000807.ge2';
onlineDir = 'C:\Users\dpqb1\Documents\Data\c103_Feb\onlineDir';
mkdir(onlineDir)
for t = 1
    b = loadGE2polar(t,fname,center,r1,r2);
    b = b./norm(b(:))*10;
%     save(fullfile(onlineDir,sprintf('polar_image_%i.mat',t)),'b')
end

% max(b(:))
% norm(b(:))
polar_image = b;

% Setup
% Parameters

P.N1 = size(polar_image,1);
P.N2 = size(polar_image,2);

% Basis function variance parameters
P.K1 = 10;
P.K2 = 15;
P.K = P.K1*P.K2;
P.sigma1 = linspace(0.5,  2,    P.K1);
P.sigma2 = linspace(0.5,  20,   P.K2);
P.basis = 'norm2';
P.mu1 = round(P.N1/2);
P.mu2 = round(P.N2/2);

% fista params
params.lambda = 0.06; % sparsity penalty
params.L = 1e3;  %
params.beta = 2; %
params.noBacktrack = 0;

params.isNonnegative = 1; % flag to enforce nonnegativity

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 100;
params.tolerance = 1e-3;

params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
params.verbose = 1;      % flag to print objective values at each iteration 
P.params = params;

% Construct dictionary
D = dictionary2D(P);
Df = fft2(D);

% Initialize solution
x_init = zeros(size(Df));
b = polar_image;

% Solve
tic
[x_hat,err,obj] = FISTA_Circulant_gpu(Df,b,params);
[awmv_rad,awmv_eta] = computeAWMV(x_hat,P);
toc


Xf = fft2(x_hat);
y_hat = ifft2(Ax_cpu(Df,x_hatf),'symmetric');

% Solve with cbpdndl for reference
% opt.AutoRho = 1;
% opt.MaxMainIter = 100;
% opt.Verbose = 1;
% opt.rho = 1;
% opt.RelaxParam = 1.8;
% opt.NonNegCoef = 1;
% X = cbpdn(D,b,0.01,opt);
% Xf = fft2(X);

y_hat = ifft2(Ax_cpu(Df,Xf),'symmetric');
figure
subplot(2,1,1)
imagesc(y_hat)
ax = gca;
ax.CLim = [0 max(b(:))];
title('Recon')

subplot(2,1,2)
imagesc(b)
ax = gca;
ax.CLim = [0 max(b(:))];
title('Data')

norm(b(:)-y_hat(:))/norm(b(:))

