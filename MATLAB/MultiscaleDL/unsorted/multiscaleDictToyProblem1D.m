% Multiscale 1D dictionary learnging toy problem
y = two_shape_problem;
y = y(1,1:25,:);

% Set up cbpdndl parameters
lambda = 0.15;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 250;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(y,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.DictFilterSizes = [1, 1, 1, 1, 1, 1;
                       25,25,25,25,25,25];
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
K = 6;
D0 = rand(1,25,K);
for i = 1:K
    D0(:,:,i) = D0(:,:,i)/norm(D0(:,:,i));
end

% Define Multiscale filter
P.N = 25;
std1 = 0.3;
a1 = 1.8;
P.stds = [std1;
          std1*a1^1;...
          std1*a1^2;...
          std1*a1^3;...
          std1*a1^4];
phi = dictionary(P);
phi(:,1) = zeros(P.N,1);
phi(1,1) = 1;
Phi = zeros([1,size(phi)]);
Phi(1,:,:) = phi;
Phif = fft2(Phi,size(y,1),size(y,2));

%% Multiscale dictionary learning

[D, X, optinf] = cbpdndl_multiscale(D0, Phi, y, lambda, opt);
Df = fft(D);
Xf = fft(X);
PhifGf = scaleDict(Phif,Df);
PhiG = ifft2(scaleDict(Phif,Df),'symmetric');
Yhat = squeeze(ifft2(sum(bsxfun(@times,PhifGf,Xf),3),'symmetric'));
plotDictionary(D)
plotMultiscaleDictionary(D,Phi)

figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat)

%% Regular dictionary learning
opt.DictFilterSizes = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1; 
                       8,8,8,8,10,10,16,16,25,25];
K = 10;
D0 = rand(1,25,K);
[D, X, optinf] = cbpdndl(D0, y, lambda, opt);
Df = fft(D);
Xf = fft(X);
Yhat = squeeze(ifft(sum(bsxfun(@times,Df,Xf),3)));
plotDictionary(D)

figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat)

function u = zpad(v, sz)
  u = zeros(sz(1), sz(2), size(v,3), size(v,4), class(v));
  u(1:size(v,1), 1:size(v,2),:,:) = v;
end