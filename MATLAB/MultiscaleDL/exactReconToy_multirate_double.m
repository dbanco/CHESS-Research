% Multiscale 1D dictionary learning toy problem
[y,AD,Dtrue,X_true,N,M,T,scales,c1,c2] = upDwn_double_multirate_problem;
y = reshape(y,[1,N,1,T]); 
% plotDataSeq(y)

[K,U] = size(scales);
V = (U-1)/2;
% plotDataSeq(y)

topDir = 'C:\Users\dpqb1\Desktop\multiDict_multirate_double\';
dName = 'doubleToy';
mkdir(topDir)

close all
% Set up cbpdndl parameters
lambda = 8e-2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.MaxCGIter = 200;
opt.CGTol = 1e-9;
opt.rho = 500*lambda + 0.5;
opt.sigma = size(y,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;


opt.DictFilterSizes = [ones(1,K);
                       M*ones(1,K)];

opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;

D0 = zeros(1,M,K);
D0(1,:,1) = AD(1,1:M,U-V); + 100*rand(1,M,1)/100;
D0(1,:,2) = AD(1,1:M,2*U-V); + 100*rand(1,M,1)/100;


%% **all below needs updating

% Check init solution (good)
opt.Y0 = zeros(size(X_true));
% opt.Y0 = X_true;
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
% 
Yhat0 = squeeze(ifft2(sum(bsxfun(@times,fft2(AD),fft2(X_true)),3),'symmetric'));
% 
f0 = figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat0)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat0,'fro')/norm(y(:),'fro')))

%% Dictionary learning
% opt.LinSolve = 'CGD';
% 
% denLim = 9;
% opt.Verbose = 0;
% [D, X, optinf, obj, relErr,output,minObj] = cbpdndlScaleSearch(Dtrue,y,lambda,U,denLim,opt);
% AD = reSampleNu(N,D,output(1),output(2),U);
% ADf = fft2(AD);
% Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));

%% Test Ax, Atb
ADtrue = reSampleNu(N,Dtrue,c1,c2,U);
plotDictUsage(ADtrue,K,1)
% ADtrue2 =  upDwnSample(Dtrue,V);
% % plotDictUsage(ADtrue(:,1:64,:)-ADtrue2,K,1)
% norm(vec(ADtrue(:,1:64,:)-ADtrue2))
% 
% Atrans1 = reSampleNuTrans2(M,X_true,c1,c2,U);
% Atrans2 = upDwnSampleTrans(X_true,V);
% imagesc(squeeze(Atrans1-Atrans2))
% norm(vec(Atrans1-Atrans2))


%% Solve
opt.LinSolve = 'CGD';
opt.Verbose = 1;
opt.MaxMainIter = 200;
[D, X, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, y, lambda, opt,c1,c2,U);
opt.numScales = 2;

AD = reSampleNu(N,D,c1,c2,U);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));

% Check against original algorithm
% [D2, X2, optinf] = cbpdndl_cg_upDwnSample(D0,y,lambda,opt);
% norm(vec(D-D2))
% norm(vec(D-Dtrue))
% norm(vec(D2-Dtrue))
% 
% norm(vec(X-X2))
% norm(vec(X-X_true))
% norm(vec(X2-X_true))
%%
% plotDictUsage(AD,K,1)
% 
f2 = figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
% saveas(f2,fullfile(topDir,['recon',dName,'.png']))
% 
f3 = figure;
subplot(2,1,1)
imagesc(squeeze(sum(sum(X,1),2)))
title('Recon')
subplot(2,1,2)
imagesc(squeeze(sum(sum(X_true,1),2)))
title('Truth')
f3.Position = [800 100 600 300];
% saveas(f3,fullfile(topDir,['dictDist',dName,'.png']))
% 
% f4 = figure;
% error_time = zeros(T,1);
% for t = 1:T
%     error_time(t) = norm(squeeze(y(:,:,:,t))'-Yhat(:,t),'fro')/norm(y(:,:,:,t),'fro');
% end
% plot(error_time)
% xlabel('time')
% ylabel('Rel Error')
% saveas(f4,fullfile(topDir,['errTime',dName,'.png']))
% 
