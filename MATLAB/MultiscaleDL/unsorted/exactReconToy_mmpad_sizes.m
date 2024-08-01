%% Multiscale 1D dictionary learning toy problem
y = loadMMPAD1D(1,1);
[N,T] = size(y);
y = reshape(y,[1,N,T]);
% plotDataSeq(y)
scales = [1;1];

U = size(scales,2);
K = 4;
M = 260;

opt = [];
% opt.DictFilterSizes = [1,1,1,1;
%                        30,M/4,M/2,M];
opt.DictFilterSizes = [1,1,1,1;
                       M,M,M,M];
% Init solution
opt.Y0 = zeros(1,N,K*U,T);
% Init dictionary
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = rand(1,M,K);
% D0(1,18:22,1:4) = 1;
% D0(1,25:35,2) = 1;
% D0(1,50:80,3) = 1;
% D0(1,20:200,4) = 1;

D0 = Pnrm(D0);

fDinit = figure;
for i = 1:K*U
    subplot(2,ceil(K*U/2),i)
    plot(D0(:,:,i),'Linewidth',1)
 set(gca, 'XtickLabel','')
set(gca, 'FontSize', 16)
end
fDinit.Position = [1 100 1800 500];
% saveas(f1,fullfile(topDir,['dict',dName,'.png']))


                   
topDir = 'C:\Users\dpqb1\Documents\Outputs\multiDict_mmpad_same_3\';
dName = sprintf('mmpad');
mkdir(topDir)
% results = struct([]);

close all
% Set up cbpdndl parameters
lambda = 3e-2;

opt.plotDict = 1;
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.MaxCGIter = 200;
opt.CGTol = 1e-9;
opt.rho = 50*lambda + 0.5;
opt.sigma = T;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
% opt.DictFilterSizes = [ones(1,K);
%                        M*ones(1,K)];

%% Dictionary learning
opt.LinSolve = 'CGD';

[D, X, optinf, obj, relErr] = cbpdndl_cg_multirate_custom(D0, y, lambda, opt, scales);
% [D, Y, optinf, Jfn, relErr] = cbpdndl_cg_multirate(D0, y, lambda, opt,2,1,3);
% save(fullfile(topDir,sprintf('output_%i.mat',i)),'D','X','opt','obj','relErr','c1','c2','output','prbCount','N','M','K','U');

% Solution
AD = reSampleCustom(N,D,scales);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));


%%

% Show dictionary
f1 = figure;
for i = 1:K*U
    subplot(2,ceil(K*U/2),i)
    plot(AD(:,:,i),'Linewidth',1)
 set(gca, 'XtickLabel','')
set(gca, 'FontSize', 16)
end
f1.Position = [1 100 1800 500];
saveas(f1,fullfile(topDir,['dict',dName,'.png']))

% Recon and data log scale
% f2 = figure;
% subplot(2,1,1)
% imagesc(log(squeeze(y)+0.05-min(y(:))))
%  set(gca, 'YtickLabel','')
% set(gca, 'FontSize', 20)
% subplot(2,1,2)
% imagesc(log(Yhat+0.05-min(Yhat(:))))
%  set(gca, 'YtickLabel','')
% set(gca, 'FontSize', 20)
%  set(gca, 'ZtickLabel','')
% title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
% saveas(f2,fullfile(topDir,['recon',dName,'.png']))

% Recon and data regular scale
f2 = figure;
subplot(2,1,1)
imagesc(squeeze(y))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
subplot(2,1,2)
imagesc(squeeze(Yhat))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
f2.Position = [1 100 800 500];
saveas(f2,fullfile(topDir,['recon',dName,'.png']))

% Recovered VDF(t)
f3 = figure;
imagesc(squeeze(sum(sum(X,1),2)))
% set(findobj(gca, 'Type', 'line'), 'LineWidth', 30)
f3.Position = [800 100 600 300];
set(gca, 'FontSize', 20)
f3.Position = [1 100 800 400];
saveas(f3,fullfile(topDir,['dictDist',dName,'.png']))

