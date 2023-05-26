center = [1025,1020];
r1 = 430;
r2 = 450;
% factor = 5;
fname = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\ff_000277.ge2';
onlineDir = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\277\';

% fname = '/nfs/chess/raw/2023-1/id1a3/miller-3528-b/c103-90-s2-4/3/ff/ff_000807.ge2';
% onlineDir = '/nfs/chess/user/dbanco/realTimeHEXRD/onlineDir';

% fname = 'C:\Users\dpqb1\Documents\Data\c103_Feb\ff_000807.ge2';
% onlineDir = 'C:\Users\dpqb1\Documents\Data\c103_Feb\onlineDir';
mkdir(onlineDir)

for t = 1:5
    b = loadGE2polar(t,fname,center,r1,r2);
    b = b./norm(b(:))*10;

    % Parameters
    P.N1 = size(b,1);
    P.N2 = size(b,2);
    
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
    params.lambda = 0.05; % sparsity penalty
    params.L = 1000;  %
    params.beta = 2; %
    params.noBacktrack = 0;
    
    params.isNonnegative = 1; % flag to enforce nonnegativity
    
    params.stoppingCriterion = 'OBJECTIVE_VALUE';
    params.maxIter = 3;
    params.tolerance = 1e-3;
    
    params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
    params.verbose = 1;      % flag to print objective values at each iteration 
    P.params = params;
    
    % Construct dictionary
    D = dictionary2D(P);
    Df = fft2(D);
    
    % Initialize solution
    x_init = zeros(P.N1,P.N2,P.K);
    % Solve
    tic
    [X,L] = FISTA_Circulant_gpu(Df,b,x_init,params);
    [awmv_rad,awmv_eta] = computeAWMV(X,P);
    x_init = X;
    Xf = fft2(X);
    b_hat = ifft2(Ax_cpu(Df,Xf),'symmetric');
    norm(b(:)-b_hat(:))/norm(b(:))
    
    toc
end

% tic
% [X,L] = FISTA_Circulant_cpu(Df,b,x_init,params);
% [awmv_rad,awmv_eta] = computeAWMV(X,P);
% 
% Xf = fft2(X);
% b_hat = ifft2(Ax_cpu(Df,Xf),'symmetric');
% norm(b(:)-b_hat(:))/norm(b(:))
% toc

% Solve with cbpdndl for reference
% opt.AutoRho = 1;
% opt.MaxMainIter = 100;
% opt.Verbose = 1;
% opt.rho = 1;
% opt.RelaxParam = 1.8;
% opt.NonNegCoef = 1;
% X = cbpdn(D,b,0.01,opt);
% Xf = fft2(X);

% y_hat = ifft2(Ax_cpu(Df,Xf),'symmetric');
% figure
% subplot(2,1,1)
% imagesc(y_hat)
% ax = gca;
% ax.CLim = [0 max(b(:))];
% title('Recon')
% 
% subplot(2,1,2)
% imagesc(b)
% ax = gca;
% ax.CLim = [0 max(b(:))];
% title('Data')



