% Multiscale 1D dictionary learnging toy problem
[y,N,T] = loadMMPAD1D;
y = reshape(y,[1,N,1,T]);

close all
% Set up cbpdndl parameters
lambda = 0.03;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.MaxCGIter = 200;
opt.CGTol = 1e-9;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(y,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
K = 1;
opt.DictFilterSizes = [ones(1,K);
                       N*ones(1,K)];
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.numScales = 5;


D0 = rand(1,N,K);
% D0(:,1:N/2,1) = 1;
% D0(1,:,2) = gaussian_basis_wrap_1D(N,N/2,5,'2-norm');
% D0(:,1:N/2,3) = 1;
% D0(1,:,4) = gaussian_basis_wrap_1D(N,N/2,5,'2-norm');
for i = 1:K
    D0(:,:,i) = D0(:,:,i)/norm(D0(:,:,i));
end

% Define Dictionary function
opt.DictSize = opt.numScales*K;
XAd = @(in,xf) ifft2(sum(bsxfun(@times,xf,fft2(decimate(in,opt.numScales))),3),'symmetric');
AtXty = @(in,xf) decimateTrans(ifft2(sum(bsxfun(@times, conj(xf), fft2(in)), 4),'symmetric'),opt.numScales);

%% Dictionary learning
opt.LinSolve = 'CGD';
[D, X, optinf] = cbpdndl_cg(D0,y,lambda,opt);
Df = fft2(D);
Xf = fft2(X);
% For regular dictionary
% Yhat = squeeze(ifft2(sum(bsxfun(@times,Df,Xf),3),'symmetric'));
% For decimated
AD = decimate(D,opt.numScales);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'));

topDir = 'C:\Users\dpqb1\Desktop\multiDict\';
dName = 'MMPAD1';
mkdir(topDir)
f1 = plotDictUsage(AD,K,opt.numScales,X);
saveas(f1,fullfile(topDir,['dictUsage_',dName,'.png']))


f2 = figure(2);
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
saveas(f2,fullfile(topDir,['recon',dName,'.png']))

f3 = figure(3);
imagesc(squeeze(sum(sum(X,1),2)))
f3.Position = [100 100 600 100];
saveas(f3,fullfile(topDir,['dictDist',dName,'.png']))

f4 = figure(4);
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