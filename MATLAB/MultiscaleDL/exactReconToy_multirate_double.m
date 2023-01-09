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
lambda = 10e-2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.MaxCGIter = 200;
opt.CGTol = 1e-8;
opt.rho = 50*lambda + 0.5;
opt.sigma = T;
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

Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = Pnrm(rand(1,M,K));
% D0(1,:,1) = AD(1,1:M,U-V); %+ 100*rand(1,M,1)/100;
% D0(1,:,2) = AD(1,1:M,2*U-V); %+ 100*rand(1,M,1)/100;


%% **all below needs updating

% Check init solution (good)
opt.Y0 = zeros(size(X_true));
% opt.Y0 = X_true;

Yhat0 = squeeze(ifft2(sum(bsxfun(@times,fft2(AD),fft2(X_true)),3),'symmetric'));
% 
f0 = figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat0)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat0,'fro')/norm(y(:),'fro')))

%% Dictionary learning
opt.LinSolve = 'CGD';

denLim = 9;
[D, X, optinf, obj, relErr,output,minObj] = cbpdndlScaleSearch(D0,y,lambda,U,denLim,opt);
opt.MaxMainIter = 200;
opt.Y0 = X;
opt.Verbose = 1;
% output = [13,4];
[D, X, optinf, obj, relErr] = cbpdndl_cg_multirate(D0,y,lambda,opt,output(1),output(2),U);

% Solution
AD = reSampleNu(N,D,output(1),output(2),U);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));

de1 = norm(squeeze(D-Dtrue),'fro')/norm(Dtrue(:),'fro');
de2 = norm(squeeze(D(:,:,[2,1])-Dtrue),'fro')/norm(Dtrue(:),'fro');
dictErr = min(de1,de2)
% results(i).obj = obj;
% results(i).relErr = norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro');
% results(i).dictErr = min(de1,de2);
% results(i).prbCount = prbCount;
% results(i).output = output;
% results(i).c1 = c1;
% results(i).c2 = c2;
% results(i).opt = opt;

ADtrue = reSampleNu(N,Dtrue,output(1),output(2),U);
f1 = figure;

for i = 1:6
    subplot(1,6,i)  
    hold off
%     plot(AD(:,:,i),'Linewidth',3)
%     hold on
    plot(ADtrue(:,:,i),'Color',[0.8500 0.3250 0.0980],'Linewidth',2)
    set(gca,'FontSize',20)
%     set(gca,'visible','off');
%     ylim([0 1])
%     if sum(squeeze(X(:,:,i,:)),'all')/sum(X(:)) > 0
%         title(sprintf('Usage: %1.2f',sum(squeeze(X(:,:,i,:)),'all')/sum(X(:)) ))
%     end
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
%     p = 1+(f.Number)*400;
    f.Position =[1 500 1000 400];
end

% plotDictUsage(AD,K,1,f1)
% saveas(f1,fullfile(topDir,['dict',dName,'.png']))

% Recon and data
f2 = figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
% saveas(f2,fullfile(topDir,['recon',dName,'.png']))

% Recovered VDF(t)
f3 = figure;
subplot(2,1,1)
imagesc(squeeze(sum(sum(X,1),2)))
title('Recon')
subplot(2,1,2)
imagesc(squeeze(sum(sum(X_true,1),2)))
title('Truth')
f3.Position = [800 100 600 300];
% saveas(f3,fullfile(topDir,['dictDist',dName,'.png']))


f4 = figure;
waterfall(squeeze(y)')
set(gca,'FontSize',20)
    set(gca,'XTickLabel',[]);
