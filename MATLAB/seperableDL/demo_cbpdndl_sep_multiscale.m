% Training images
close all
clear all
y = loadMMPAD_eta_om();

%% Construct initial dictionary
[N1,N2,T] = size(y);
K=2; M1 = 8; M2 = 12;
D0v = zeros(M1,1,K);
D0v(:,:,1) = gaussian_basis_wrap_1D(M1,M1/2,2,'2-norm');
D0v(:,:,2) = gaussian_basis_wrap_1D(M1,M1/2,2,'2-norm');

D0h = zeros(1,M2,K);
D0h(:,:,1) = gaussian_basis_wrap_1D(M2,M2/2,6,'2-norm');
D0h(:,:,2) = gaussian_basis_wrap_1D(M2,M2/2,6,'2-norm');
  
for c=1:K
   D0(:,:,c) = D0v(:,:,c)*D0h(:,:,c);
end

% Set up cbpdndl parameters
lambda = 8e-2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 4;
opt.MaxCGIter = 200;

opt.rho = (500*lambda + 0.5);
opt.sigma = (size(y,3));
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.LinSolve = 'CG';
opt.CGTol = 1e-3;
opt.CGTolAuto = 0;
opt.NonnegativeDict = 1;
opt.NonNegCoef = 1;

% Do dictionary learning
c1 = 2;
c2 = 1;
Ufactors = 3;
[D,AG, Y, Gh, Gv, optinf, Jfn, relErr] = cbpdndl_cg_multirate_sep(D0,D0h,D0v, y, lambda, opt, c1,c2,Ufactors);

AGh = reSampleNu2d(N2,Gh,c1,c2,Ufactors);
AGv = reSampleNu2d(N1,Gv,c1,c2,Ufactors);
[AG,NormVals] = combineDict(N1,N2,AGh,AGv,K,Ufactors);

yhat = squeeze(ifft2(sum(bsxfun(@times,fft2(AG),fft2(Y)),3),'symmetric'));

figure;
imdisp(tiledict(AG));

figure;
subplot(2,2,1)
imagesc(y(:,:,1))
subplot(2,2,2)
imagesc(y(:,:,20))
subplot(2,2,3)
imagesc(yhat(:,:,1))
subplot(2,2,4)
imagesc(yhat(:,:,20))

% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');
