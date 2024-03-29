%% Multiscale 1D dictionary learning toy problem
[y,N,M,T] = gaussian_pwlinear_2to10_problem;
y = reshape(y,[1,N,T]);

topDir = 'C:\Users\dpqb1\Documents\Outputs\multiScales_toy2_many1\';
dName = sprintf('gs_sq');
mkdir(topDir)

plotDataSeq(y,topDir)
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],8,100, 1/8);
scales{2} = genRationals([0;1],[1;1],8,100, 1/8);
% scales{3} = genRationals([0;1],[1;1],8,100, 1/8);

Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
Utotal = sum(Uarray);


opt = [];
% opt.DictFilterSizes = [1,1,1;
%                        M,M,M];
opt.DictFilterSizes = [1,1;
                       M,M];
% Init solution
opt.Y0 = zeros(1,N,Utotal,T);
% Init dictionary
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = zeros(1,M,K);
D0(1,60:140,1) = 1;
D0(1,40:180,2) = 1;
% D0(1,20:200,3) = 1;
D0 = Pnrm(D0);
                  
close all
% Set up cbpdndl parameters
lambda = 7e-2;

opt.plotDict = 0;
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

%% Dictionary learning
opt.LinSolve = 'CGD';

[D, X, optinf, obj, relErr] = cbpdndl_cg_multiScales(D0, y, lambda, opt, scales);
% save(fullfile(topDir,sprintf('output_%i.mat',i)),'D','X','opt','obj','relErr','c1','c2','output','prbCount','N','M','K','U');

% Solution
AD = reSampleCustomArray(N,D,scales);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));


%%

% Show dictionary
f1 = figure;
for i = 1:Utotal
    subplot(7,7,i)
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
saveas(f3,fullfile(topDir,['vdf',dName,'.png']))

% Spatial placement of atom groups
i = 1;
for k = 1:K
    figure;
    Ui = size(scales{k},2) + i - 1;
    imagesc( squeeze(sum(X(:,:,i:Ui,:),3)) )
    i = i + Ui;
end

% % Spatial placement of atom groups
% for k = 1:Utotal
%     xk = squeeze(X(:,:,k,:));
%     if(max(xk(:))>1e-8)
%         figure(k);
%         imagesc( xk )
%     end
% end
