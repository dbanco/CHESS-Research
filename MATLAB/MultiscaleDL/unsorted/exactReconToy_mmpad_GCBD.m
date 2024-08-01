%% Multiscale 1D dictionary learning toy problem
y = loadMMPAD1D(1,1);
[N,T] = size(y);
y = reshape(y,[1,N,T]);
% plotDataSeq(y)
% scales = [[1;8],[1;6],[1;5],[1;4],[1;3],[1;2],[2;3],[1;1],[3;2],...
%           [2;1],[3;1],[4;1],[5;1],[6;1],[7;1],[8;1],[49;6]];
% scales = [[1;20],[1;18],[1;16],[1;10],[1;8],[1;6],[1;5],[1;4],[1;3],[1;2],...
%           [5;9],[4;7],[3;5],[5;8],[2;3],[5;7],[3;4],[4;5],[1;1]];
% scales = [[1;2],[1;1],[2;1]];
scales = cell(3,1);
scales{1} = genRationals([1;1],[2;1],3,100, 1/1);
scales{2} = genRationals([0;1],[2;1],3,100, 1/2);
scales{3} = genRationals([0;1],[1;1],3,100, 1/3);

Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
Utotal = sum(Uarray);
K = 6;
M = 50;

opt = [];
opt.DictFilterSizes = [1,1,1,1,1,1;
                       6,12,18,24,36,50];
% Init solution
opt.Y0 = zeros(1,N,K,T);
% Init dictionary
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = zeros(1,M,K);
D0(1,1:6,1) = 1;
D0(1,1:12,2) = 1;
D0(1,1:18,3) = 1;
D0(1,1:24,4) = 1;
D0(1,1:36,5) = 1;
D0(1,1:50,6) = 1;
D0 = Pnrm(D0);
                   
topDir = 'C:\Users\dpqb1\Documents\Outputs\multiDict_GCBD\';
dName = sprintf('mmpad');
mkdir(topDir)
% results = struct([]);

close all
% Set up cbpdndl parameters
lambda = 7e-2;

% opt.showImage = 1;
opt.epsilon = 5;
opt.delta = 1e-6;
opt.Verbose = 1;
opt.MaxMainIter = 50;
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

% Normalisation projection
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

% Projection of filter to full image size and its transpose
% (zero-pad and crop respectively)
Pzp = @(x) zpad(x, [1,N]);
PzpT = @(x) bcrop(x, opt.DictFilterSizes);
Pnonneg = @(x) makeNonneg(x);

% Projection of dictionary filters onto constraint set
if opt.NonnegativeDict
  Pcn = @(x) removeNan(Pnrm(Pzp(PzpT(Pnonneg(x)))));
else
  Pcn = @(x) Pnrm(Pzp(PzpT(x)));
end

% X0 = zeros(1,N,3,T);
% X0(1,1,1,1) = 1;
% X0f = fft2(X0);
% Df = fft2(Pcn(D0));
% Dx = ifft2(sum(bsxfun(@times,Df,X0f),3),'symmetric');
% plot(Dx(1,:,1,1))

%% Dictionary learning
opt.LinSolve = 'CGD';

[D, X, optinf] = cbpdndl_GCBD(D0, y, lambda, opt);
% [D, Y, optinf, Jfn, relErr] = cbpdndl_cg_multirate(D0, y, lambda, opt,2,1,3);
% save(fullfile(topDir,sprintf('output_%i.mat',i)),'D','X','opt','obj','relErr','c1','c2','output','prbCount','N','M','K','U');

% Solution
% AD = reSampleCustomArray(N,D,scales);
Df = fft2(Pzp(D));
Yhat = squeeze(ifft2(sum(bsxfun(@times,Df,fft2(X)),3),'symmetric'));


%%

% Show dictionary
f1 = figure;
for i = 1:6
    subplot(3,2,i)
    plot(D(:,:,i),'Linewidth',1)
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

function y = makeNonneg(x)
    y=x;
    y(x<0)=0;
end

function y = removeNan(x)
    y=x;
    y(isnan(x))=0;
end

function u = zpad(v, sz)

  u = zeros(sz(1), sz(2), size(v,3), size(v,4), class(v));
  u(1:size(v,1), 1:size(v,2),:,:) = v;

end

function u = bcrop(v, sz)

  if numel(sz) <= 2,
    if numel(sz) == 1
      cs = [sz sz];
    else
      cs = sz;
    end
    u = v(1:cs(1), 1:cs(2), :);
  else
    cs = max(sz,[],2);
    u = zeros(cs(1), cs(2), size(v,3), class(v));
    for k = 1:size(v,3),
      u(1:sz(1,k), 1:sz(2,k), k) = v(1:sz(1,k), 1:sz(2,k), k);
    end
  end

end