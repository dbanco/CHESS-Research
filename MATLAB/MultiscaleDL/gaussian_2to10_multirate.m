% Multiscale 1D dictionary learning toy problem
[y,N,M,T] = gaussian_2to10_problem;
M = 32;
y = reshape(y,[1,N,1,T]); 
K = 1;
U = 5;
V = (U-1)/2;
% plotDataSeq(y)

close all
% Set up cbpdndl parameters
lambda = 8e-2;
opt = [];
opt.Verbose = 1;
opt.plotDict = 0;
opt.MaxMainIter = 200;
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

D0 = zeros(1,M,K);
for k = 1:K
    D0(1,round(M/3):round(2*M/3),k) = rand;
    D0(:,:,k) = D0(:,:,k)/norm(D0(:,:,k));
end

% D0(1,:,2) = AD(1,1:M,2*U-V); + 100*rand(1,M,1)/100;

%% Solve
opt.Y0 = zeros(1,N,K*U,T);
opt.LinSolve = 'CGD';
opt.Verbose = 1;
denLim = 11;

[D, X, optinf, obj, relErr,output,minObj] = cbpdndlScaleSearch(D0,y,lambda,U,denLim,opt);
c1 = output(1);
c2 = output(2);

% c1 = 13;
% c2 = 7;
% [D, X, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, y, lambda, opt,c1,c2,U);
 

AD = reSampleNu(N,D,c1,c2,U);

%% Test norms
% c1c2 = [2:2:64,231;
%         1:32,111];
% [n1,n2] = size(c1c2);
% normLo = zeros(16,1);
% normHi = zeros(16,1);
% Dtest = randn(1,32,1);
% for i = 1:n2
%     AD = reSampleNu(N,Dtest,c1c2(1,i),c1c2(2,i),U);
%     normLo(i) = norm(AD(1,:,1));
%     normHi(i) = norm(AD(1,:,3));
% end
% figure
% hold on
% plot(normLo,'o-')
% plot(normHi,'o-')
% x = c1c2(1,:);
% plot( 1.5./(1+exp(-0.7*x)),'o-')

%%
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));

topDir = 'C:\Users\dpqb1\Desktop\gaussian_2to10_problem_K2\';
dName = 'gaussian_2to10';
mkdir(topDir)

% Show dicitonary
[~,~,KU] = size(AD);
for i = 1:KU
    subplot(K,KU/K,i)
    plot(AD(:,:,i),'Linewidth',2)
if i > 3    
    ylim([0 0.52])
else
    ylim([0 0.52])
end
    set(gca,'XTickLabel',[]);
end

% saveas(f1,fullfile(topDir,['dict',dName,'.png']))

% Show usage
f2 = figure;
imagesc(squeeze(sum(sum(X,1),2)))
% title('vdf')
set(gca, 'FontSize', 20)
saveas(f2,fullfile(topDir,['vdf',dName,'.png']))

% Show recon
f3 = figure;
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
saveas(f3,fullfile(topDir,['recon',dName,'.png']))

f4 = figure;
waterfall(squeeze(y(1,80:228,1,:))')
set(gca, 'XtickLabel','')
set(gca, 'ZtickLabel','')
set(gca, 'FontSize', 24)