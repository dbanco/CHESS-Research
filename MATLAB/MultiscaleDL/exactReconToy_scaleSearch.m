% Multiscale 1D dictionary learnging toy problem
[y,yd,X_true,N,T] = upDwn_double_problem;
y = reshape(y,[1,N,1,T]); 
% plotDataSeq(y)

topDir = 'C:\Users\dpqb1\Desktop\multiDict_upDwn_double\';
dName = 'doubleToy';
mkdir(topDir)

U = 2;
close all
% Set up cbpdndl parameters
lambda = 3e-1;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 500;
opt.MaxCGIter = 200;
opt.CGTol = 1e-9;
opt.rho = 70*lambda + 0.5;
opt.sigma = size(y,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
K = 2;
N0 = N/2^U;
opt.DictFilterSizes = [ones(1,K);
                       N0*ones(1,K)];
%                        N,N,N,N,N,N];
%                        N,N/2,N/4,N,N/2,N/4];


opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.numScales = U;

D0 = zeros(1,N0,K);
% D0(:,1:N/2,1) = 1;
% D0(1,:,2) = gaussian_basis_wrap_1D(N,N/2,5,'2-norm');
% D0(:,1:N/2,3) = 1;
% D0(1,:,4) = gaussian_basis_wrap_1D(N,N/2,5,'2-norm');
% for i = 1:K
%     D0(:,:,i) = D0(:,:,i)/norm(D0(:,:,i));
% end
D0(1,:,1) = yd(1,1:N0,U+1) + 100*rand(1,N0,1)/100;
D0(1,:,2) = yd(1,1:N0,U+1+2*U+1) + 100*rand(1,N0,1)/100;

% Initial Solution
fff=checkDict_upDwnSample(X_true,D0,opt.numScales);
saveas(fff,fullfile(topDir,['init',dName,'.png']))

% Show true dictionary
f5 = checkDict_upDwnSample(X_true,yd(1,1:N0,[U+1,U+1+2*U+1 ]),opt.numScales);
saveas(f5,fullfile(topDir,['truth',dName,'.png']))

% Define Dictionary function
opt.DictSize = (2*U+1)*K;
XAd = @(in,xf) ifft2(sum(bsxfun(@times,xf,fft2(upDwnSample(in,opt.numScales))),3),'symmetric');
AtXty = @(in,xf) upDwnSampleTrans(ifft2(sum(bsxfun(@times, conj(xf), fft2(in)), 4),'symmetric'),opt.numScales);

% Check init solution (good)
opt.Y0 = zeros(size(X_true));
% opt.Y0 = X_true;
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
% 
AD0 = upDwnSample(D0,opt.numScales);
AD0f = fft2(AD0);
Yhat0 = squeeze(ifft2(sum(bsxfun(@times,AD0f,fft2(opt.Y0)),3),'symmetric'));

f0 = figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat0)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat0,'fro')/norm(y(:),'fro')))

%% Dictionary learning
opt.LinSolve = 'CGD';
[D, X, optinf] = cbpdndlScaleSearch(y,N,D0,opt.Y0,U,denLim,opt);
Df = fft2(D);
Xf = fft2(X);
% For regular dictionary
% Yhat = squeeze(ifft2(sum(bsxfun(@times,Df,Xf),3),'symmetric'));
% For decimated
AD = upDwnSample(D,opt.numScales);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'));

f1 = checkDict_upDwnSample(X,D,opt.numScales);
% plotDictionary(decimate(D0,opt.numScales))
saveas(f1,fullfile(topDir,['dictUsage_',dName,'.png']))

f2 = figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
saveas(f2,fullfile(topDir,['recon',dName,'.png']))

f3 = figure;
subplot(2,1,1)
imagesc(squeeze(sum(sum(X,1),2)))
title('Recon')
subplot(2,1,2)
imagesc(squeeze(sum(sum(X_true,1),2)))
title('Truth')
f3.Position = [800 100 600 300];
saveas(f3,fullfile(topDir,['dictDist',dName,'.png']))

f4 = figure;
error_time = zeros(T,1);
for t = 1:T
    error_time(t) = norm(squeeze(y(:,:,:,t))'-Yhat(:,t),'fro')/norm(y(:,:,:,t),'fro');
end
plot(error_time)
xlabel('time')
ylabel('Rel Error')
saveas(f4,fullfile(topDir,['errTime',dName,'.png']))

% function u = zpad(v, sz)
%   u = zeros(sz(1), sz(2), size(v,3), size(v,4), class(v));
%   u(1:size(v,1), 1:size(v,2),:,:) = v;
% end