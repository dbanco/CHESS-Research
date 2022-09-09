% Multiscale 1D dictionary learnging toy problem
y = two_shape_problem;
close all
% Set up cbpdndl parameters
lambda = 0.01;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(y,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.DictFilterSizes = [1, 1;
                       64,64];
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.numScales = 1;

K = 2;
D0 = zeros(1,64,K);
D0(:,1:32,1) = 1;
D0(1,:,2) = gaussian_basis_wrap_1D(64,32,5,'2-norm');
for i = 1:K
    D0(:,:,i) = D0(:,:,i)/norm(D0(:,:,i));
end

% Define Dictionary function
opt.DictSize = opt.numScales*K;
XAd = @(in,xf) ifft2(sum(bsxfun(@times,xf,fft2(decimate(in,opt.numScales))),3),'symmetric');
AtXty = @(in,xf) decimateTrans(ifft2(sum(bsxfun(@times, xf, fft2(in)), 4),'symmetric'),opt.numScales);
%% Multiscale dictionary learning
[D, X, optinf] = cbpdndl_functional(D0,y,lambda,opt);
Xf = fft2(X);
AD = decimate(D,opt.numScales);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'));
AD0 = decimate(D0,opt.numScales);

figure(1)
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat)
plotDictionary(AD,K,opt.numScales);
plotDictionary(AD0,K,opt.numScales);

% %% Regular dictionary learning
% opt.DictFilterSizes = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1; 
%                        8,8,8,8,10,10,16,16,25,25];
% K = 10;
% D0 = rand(1,25,K);
% [D, X, optinf] = cbpdndl(D0, y, lambda, opt);
% Df = fft(D);
% Xf = fft(X);
% Yhat = squeeze(ifft(sum(bsxfun(@times,Df,Xf),3)));
% plotDictionary(D)
% 
% figure;
% subplot(1,2,1)
% imagesc(squeeze(y))
% subplot(1,2,2)
% imagesc(Yhat)

function u = zpad(v, sz)
  u = zeros(sz(1), sz(2), size(v,3), size(v,4), class(v));
  u(1:size(v,1), 1:size(v,2),:,:) = v;
end