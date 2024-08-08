% Multiscale 1D dictionary learning mmpad problem
[y,N,T] = loadMMPAD1D;
y = reshape(y,[1,N,1,T]);
plotDataSeq(y)
%%
M = 32;
K = 2;
U = 3;
V = (U-1)/2;
c1s = [2 3 3 5 7];
c2s = [1 1 2 3 2];

Dout = cell(5,1);
Xout = cell(5,1);

for i = 1:5
    c1 = c1s(i);
    c2 = c2s(i);


    topDir = 'C:\Users\dpqb1\Desktop\multiDict_multirate_double\';
    dName = 'doubleToy';
    mkdir(topDir)

    close all
    % Set up cbpdndl parameters
    lambda = 8e-2;
    opt = [];
    opt.Verbose = 1;
    opt.MaxMainIter = 10;
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

    opt.DictFilterSizes = [ones(1,K);
                           M*ones(1,K)];

    opt.NonNegCoef = 1;
    opt.NonnegativeDict = 1;

    D0 = zeros(1,M,K);
    D0(1,:,1) = rand(1,M,1);
    D0(1,:,2) = rand(1,M,1);

    opt.Y0 = zeros(1,N,K*U,T);

    % Yhat0 = squeeze(ifft2(sum(bsxfun(@times,fft2(AD),fft2(X_true)),3),'symmetric')); 
    % f0 = figure;
    % subplot(1,2,1)
    % imagesc(squeeze(y))
    % subplot(1,2,2)
    % imagesc(Yhat0)
    % title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat0,'fro')/norm(y(:),'fro')))
    
% Solve
    opt.LinSolve = 'CGD';
    opt.Verbose = 1;
    opt.MaxMainIter = 200;
    [D, X, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, y, lambda, opt,c1,c2,U);
    Dout{i} = D;
    Xout{i} = X;
end

%%
for i = 1:5
    AD = reSampleNu(N,Dout{i},c1s(i),c2s(i),U);
    ADf = fft2(AD);
    Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(Xout{i})),3),'symmetric'));
    
    % Show dicitonary
    plotDictUsage(AD,K,1)

    % Show usage
    f2 = figure;
    imagesc(squeeze(sum(sum(X,1),2)))
    title('Recon')

    % Show recon
    f3 = figure;
    subplot(1,2,1)
    imagesc(squeeze(y))
    subplot(1,2,2)
    imagesc(Yhat)
    title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
    pause()
end
%     saveas(f2,fullfile(topDir,['recon',dName,'.png']))
 
% f3.Position = [800 100 600 300];
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
