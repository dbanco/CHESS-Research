lambdaVals = logspace(-2,0,150);
lambdaRegVals = [0 logspace(-2,2,99)];

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 300;
opt.MaxCGIter = 100;
opt.CGTol = 1e-6;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
% Rho1 and Rho2 params  
opt.rho1 = 50;
opt.rho2 = 5;
opt.AutoRho1 = 1;
opt.AutoRho1Period = 1;
opt.AutoRho2 = 1;
opt.AutoRho2Period = 1;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.a_min = 1e-4;
opt.lambda_min = 1e-4;
opt.adapt_a = false;
opt.adapt_lambda = false;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.useGpu = false;
opt.Xfixed = 0;
opt.Dfixed = 0;
opt.Recenter = 0;
opt.a = 1;
opt.useMin = false;
opt.AdaptIters = 100;
opt.a_via_lam = true;
opt.l1_iters = 10;
opt.mcdl_init = false;
opt.ism_init = false;
opt.L = 1;

% Multiscale dictionary setup
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);

funcName = 'sim_mcdl_reg_wrapper';
k = 1;

datasets = {'pseudo-voigt_unmatched'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenter = {0,1};

% Noise level
dataset = datasets{1};
sig_ind = 1:5;
SNRs = [20,16,12,8,4];
sigmas = zeros(numel(SNRs),1);
for i = sig_ind
    sigmas(i) = SNRtoSigma(SNRs(i),dataset);
end

ind1 = 1:numel(lambdaVals);

% Regularizer
opt.regularizer = 'softmin';
a_n = 1;


% ISM vs Gradient Descent comparison

% Data  
rng(1);
[y,~,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigmas(4),dataset);
y = reshape(y,[1,N,1,T]);


% Model Setup
J = size(scales{1},2);
opt.DictFilterSizes = [1; M];
center = (M+1)/2;

lambda = lambdaVals(90);
lambda2 = lambdaRegVals(1);

opt.lambda = lambda;
opt.lambda2 = lambda2;

% Initialize Coefficients and Dictionary
opt.coefInit = 'zeros';
opt.dictInit = 'flat';
opt = initXD(opt,N,M,K,J,T,Xtrue,Dtrue);

Spad = padarray(y,[0 M-1 0 0],0,'pre');
Spadf = fft2(Spad);
[AG,NormVals,Shifts] = reSampleCustomArrayCenter3(N,Dtrue,scales,center);
AG = padarray(AG,[0 M-1 0 0],0,'post');
AGf = fft2(AG);
Aop = @(xf) sum(bsxfun(@times,AGf,xf),3);
Atop = @(rf) bsxfun(@times, conj(AGf), rf);

AGSf = Atop(Spadf);

Y = zeros(size(Xtrue));
U = Y;
U1 = U;
U2 = U;
U3 = U;
Yf = fft2(Y);

%% Iteration 1
bf = AGSf + opt.rho1*fft2(-Y+U);
opts.b = 0;
opts.L = 1;

Xf2 = solvedbi_sm(AGf, opt.rho1, bf);
X2 = ifft2(Xf2,'symmetric');

options = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'HessUpdate','lbfgs', ...  
    'SpecifyObjectiveGradient',true, ...
    'MaxIterations',200, ...
    'Display','iter', ...
    'OptimalityTolerance',1e-16, ...
    'StepTolerance',1e-16);


x0 = Y(:);
[xopt, fval, exitflag, output] = fminunc(@(x)obj_and_grad_x(x,AGf,Spadf,Y,U,opt.rho1,opt.lambda2,K,J), x0, options);
X3 = reshape(xopt, size(Y));
Xf3 = fft2(X3);

NUM_ITERS = 6;

recon2 = squeeze(unpad(ifft2(sum(bsxfun(@times,AGf,Xf2),3),'symmetric'),M-1,'pre'));
recon3 = squeeze(unpad(ifft2(sum(bsxfun(@times,AGf,Xf3),3),'symmetric'),M-1,'pre'));

subplot(3,NUM_ITERS,1)
hold on
plot(recon2(:,1))
plot(recon3(:,1))
title('Iter 1, Time 1')
legend('ISM','LBFGS')

subplot(3,NUM_ITERS,NUM_ITERS+1)
hold on
plot(recon2(:,15))
plot(recon3(:,15))
title('Iter 1, Time 15')

subplot(3,NUM_ITERS,2*NUM_ITERS+1)
hold on
plot(recon2(:,30))
plot(recon3(:,30))
title('Iter 1, Time 30')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % flatten helpers
% 
% Xf0 = Yf;
% vec = @(x) x(:);
% reshape_x = @(v) reshape(v, size(Xf0));  % Xf0 is initial freq coeffs
% 
% % operator: u_vec -> (Atop(Aop(reshape(u_vec))) + rho * reshape(ifft2(u_vec)) fft?) 
% % More efficient if operate in freq domain mostly:
% AhAvop = @(u_vec) vec( bsxfun(@times, conj(AGf), sum(bsxfun(@times, AGf, reshape(u_vec, size(Xf0))),3) ) ) + opt.rho1 * u_vec;
% 
% % Preconditioner (diagonal, positive)
% Pdiag = sum(abs(AGf).^2, 3) + opt.rho1;  % [Nx,Ny]
% Pdiag = repmat(Pdiag, [1,1,J,T]);     % broadcast to all channels
% Mfun  = @(v) v ./ vec(Pdiag);
% 
% 
% v = randn(size(bvec));
% y = AhAvop(v);
% z = Mfun(v);
% assert(isequal(size(y), size(bvec)));
% assert(isequal(size(z), size(bvec)));
% 
% bvec = vec(bf);           % frequency RHS flattened
% tol = 1e-6; maxit = 500;
% [xv,flag,relres,iter] = pcg(@(u)AhAvop(u), bf(:), tol, maxit, @(v)Mfun(v), [], vec(Xf0));
% Xf4 = reshape(xv, size(Xf0));
% X4 = ifft2(Xf4,'symmetric');
% 
% Pdiag2 = sum(abs(AGf).^2, 3) + opt.rho1;
% recon4 = squeeze(unpad(ifft2(sum(bsxfun(@times,AGf,Xf4),3),'symmetric'),M-1,'pre'));
% subplot(3,3,1)
% hold on
% plot(recon4(:,1))
% title('Iter 1, Time 1')
% legend('GD','ISM','LBFGS','PCG')
% 
% subplot(3,3,4)
% hold on
% plot(recon4(:,15))
% title('Iter 1, Time 15')
% 
% subplot(3,3,7)
% hold on
% plot(recon4(:,30))
% title('Iter 1, Time 30')


%% Iteration 2

for iter = 2:5

    opt.Penalty = 'l1-norm';
    opt.L1Weight = 1;
    opt.NoBndryCross = false;
    % Y1 = apply_sparse_penalty(X1,opt.lambda,k,opt,a_n);
    Y2 = apply_sparse_penalty(X2,opt.lambda,k,opt,a_n);
    Y3 = apply_sparse_penalty(X3,opt.lambda,k,opt,a_n);
    
    Yf2 = fft2(Y2);
    Yf3 = fft2(Y3);
    
    disp(norm(Yf2(:)-Yf3(:)))

    U2 = U2 + X2 - Y2;
    U3 = U3 + X3 - Y3;
    
    bf2 = AGSf + opt.rho1*fft2(Y2-U2);
    bf3 = AGSf + opt.rho1*fft2(Y3-U3);
    
    Xf2 = solvedbi_sm(AGf, opt.rho1, bf2);
    [xopt, fval, exitflag, output] = fminunc(@(x)obj_and_grad_x(x,AGf,Spadf,Y3,U3,opt.rho1,opt.lambda2,K,J), Y3(:), options);
    X_opt = reshape(xopt, size(Y));
    Xf3 = fft2(X_opt);
    
    recon2 = squeeze(unpad(ifft2(sum(bsxfun(@times,AGf,Xf2),3),'symmetric'),M-1,'pre'));
    recon3 = squeeze(unpad(ifft2(sum(bsxfun(@times,AGf,Xf3),3),'symmetric'),M-1,'pre'));
    
    subplot(3,NUM_ITERS,iter)
    plot(recon2(:,1))
    hold on
    plot(recon3(:,1))
    title('Iter 2, Time 1')
    legend('ISM','LBFGS')
    
    subplot(3,NUM_ITERS,NUM_ITERS+iter)
    plot(recon2(:,15))
    hold on
    plot(recon3(:,15))
    title('Iter 2, Time 15')
    
    subplot(3,NUM_ITERS,2*NUM_ITERS+iter)
    plot(recon2(:,30))
    hold on
    plot(recon3(:,30))
    title('Iter 2, Time 30')
end