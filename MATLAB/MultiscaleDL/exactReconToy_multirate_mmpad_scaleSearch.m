% Multiscale 1D dictionary learning mmpad problem

[y,N,T] = loadMMPAD1D;
y = reshape(y,[1,N,1,T]);
% plotDataSeq(y)
M = 32;
results = struct([]);
for K = 2:5
    U = 3;
    V = (U-1)/2;
    c1s = [2 3 3 5 7];
    c2s = [1 1 2 3 2];

    Dout = cell(5,1);
    Xout = cell(5,1);

    close all
    % Set up cbpdndl parameters
    lambda = 8e-3;
    opt = [];
    opt.Verbose = 1;
    opt.MaxMainIter = 40;
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
    for k = 1:K
        D0(1,:,k) = rand(1,M,1);
    end

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
    denLim = 7;
%     [D, Y, optinf, obj, relErr,output,minObj] = cbpdndlScaleSearch(D0,y,lambda,U,denLim,opt);
[D, Y, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, y, lambda, opt,4,1,U);
        
    results(K-1).D = D;
    results(K-1).Y = Y;
    results(K-1).obj = obj;
%     results(K-1).output = output;
%     results(K-1).minObj = minObj;
end

%%  

topDir = 'C:\Users\dpqb1\Desktop\scaleSearch_K_mmpad128_Norm_4_1\';
mkdir(topDir)
for K = 2:5
    close all
    dName = sprintf('_K_%i',K);
    i = K-1;
%     c1 = results(i).output(1);
%     c2 = results(i).output(2);
    D = results(i).D;
    X = results(i).Y;

    AD = reSampleNu(N,D,c1,c2,U);
    ADf = fft2(AD);
    Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));
    
    % Show dicitonary
    f1 = plotDictUsage(AD,K,1);
    saveas(f1,fullfile(topDir,['dict',dName,'.png']))

    % Show usage
    f2 = figure;
    imagesc(squeeze(sum(sum(X,1),2)))
    title('vdf')
    saveas(f2,fullfile(topDir,['vdf',dName,'.png']))
    
    % Show recon
    f3 = figure;
    subplot(1,2,1)
    imagesc(squeeze(y))
    subplot(1,2,2)
    imagesc(Yhat)
    title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
    saveas(f3,fullfile(topDir,['recon',dName,'.png']))
end

 
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
