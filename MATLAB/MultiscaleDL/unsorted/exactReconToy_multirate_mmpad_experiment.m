%% Multiscale 1D dictionary learning toy problem
y = loadMMPAD1D(1,1);
[N,T] = size(y);
y = reshape(y,[1,N,T]);
% plotDataSeq(y)
U = 3;
K = 2;
M = 64;
opt.DictFilterSizes = [1 1;
                      32 64 ];
lowLim = [1;1];
hiLim = [4;1]; 
denLim = 9;
V = (U-1)/2;

% topDir = 'C:\Users\dpqb1\Desktop\multiDict_multirate_mmpad_experiment\';
% dName = sprintf('mmpad_%i_c1c2_%i_%i',i,c1,c2);
% mkdir(topDir)
% results = struct([]);

close all
% Set up cbpdndl parameters
lambda = 3e-2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 30;
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

Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = Pnrm(rand(1,M,K));
% D0(1,:,1) = AD(1,1:M,U-V); %+ 100*rand(1,M,1)/100;
% D0(1,:,2) = AD(1,1:M,2*U-V); %+ 100*rand(1,M,1)/100;

% Init solution 
opt.Y0 = zeros(1,N,K*U,T);
% opt.Y0 = X_true;

% View init solution
% Yhat0 = squeeze(ifft2(sum(bsxfun(@times,fft2(AD),fft2(X_true)),3),'symmetric'));
% f0 = figure;
% subplot(1,2,1)
% imagesc(squeeze(y))
% subplot(1,2,2)
% imagesc(Yhat0)
% title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat0,'fro')/norm(y(:),'fro')))

%% Dictionary learning
opt.LinSolve = 'CGD';
% [D, X, optinf, obj, relErr,output,minObj,prbCount] = cbpdndlScaleSearch(D0,y,lambda,U,denLim,opt);
opt.MaxMainIter = 200;
output = [3;1];
[D, X, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, y, lambda, opt,output(1),output(2),U);
% save(fullfile(topDir,sprintf('output_%i.mat',i)),'D','X','opt','obj','relErr','c1','c2','output','prbCount','N','M','K','U');

% Solution
AD = reSampleNu(N,D,output(1),output(2),U);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));


% results(i).obj = obj;
% results(i).relErr = norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro');
% de1 = norm(squeeze(D-Dtrue),'fro')/norm(Dtrue(:),'fro');
% de2 = norm(squeeze(D(:,:,[2,1])-Dtrue),'fro')/norm(Dtrue(:),'fro');
% results(i).dictErr = min(de1,de2);
% results(i).prbCount = prbCount;
% results(i).output = output;
% results(i).c1 = c1;
% results(i).c2 = c2;
% results(i).opt = opt;

%%
f1 = figure;
for i = 1:6
    subplot(1,6,i)
    plot(AD(:,:,i),'Linewidth',2)
 set(gca, 'XtickLabel','')
set(gca, 'FontSize', 16)
end

% saveas(f1,fullfile(topDir,['dict',dName,'.png']))

% Recon and data
f2 = figure;
subplot(2,1,1)
imagesc(log(squeeze(y)+0.05-min(y(:))))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
subplot(2,1,2)
imagesc(log(Yhat+0.05-min(Yhat(:))))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
%  set(gca, 'ZtickLabel','')
% title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
% saveas(f2,fullfile(topDir,['recon',dName,'.png']))

f22 = figure;
subplot(2,1,1)
imagesc(squeeze(y))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
subplot(2,1,2)
imagesc(squeeze(Yhat))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
% title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
% saveas(f2,fullfile(topDir,['recon',dName,'.png']))

f23 = figure;
subplot(2,1,1)
imagesc(squeeze(y))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
subplot(2,1,2)
imagesc(squeeze(Yhat))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
% title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
% saveas(f2,fullfile(topDir,['recon',dName,'.png']))

% Recovered VDF(t)
f3 = figure;
waterfall(squeeze(sum(sum(X,1),2)), 'LineWidth', 30)
% set(findobj(gca, 'Type', 'line'), 'LineWidth', 30)
f3.Position = [800 100 600 300];
set(gca, 'FontSize', 20)
% saveas(f3,fullfile(topDir,['dictDist',dName,'.png']))

% Recon and data
f4 = figure;
imagesc(log(squeeze(y)+0.05-min(y(:))))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)


% save(fullfile(topDir,'results_exp2.mat'),'results')
