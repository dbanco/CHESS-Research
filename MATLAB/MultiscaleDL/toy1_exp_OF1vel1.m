%% Multiscale 1D dictionary learning toy problem
% Directory
% topDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1';
topDir = '/cluster/home/dbanco02/Outputs/toy1_exp_OF1vel1';
% mkdir(topDir)

% Experiment Setup
% lambdas = logspace(-2,-0.5,20);
lambdas = [logspace(-2,-0.5,20), logspace(-0.4,0,10)];
sigmas = 0.01:0.01:0.05;

% Data parameters
[y,~,N,M,T] = gaus_linear_osc_signal(0.04);
y = reshape(y,[1,N,T]);

% Model Setup
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],8,100, 1/8);
scales{2} = genRationals([0;1],[1;1],8,100, 1/8);
Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
KM = sum(Uarray);
opt = [];
opt.DictFilterSizes = [1,1;...
                       M,M];

% Init solution
opt.Y0 = zeros(1,N,KM,T);
opt.Z0 = zeros(1,1,KM,T-1);

% Init dictionary
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = zeros(1,M,K);
D0(1,121:180,1) = 1;
D0(1,101:200,2) = 1;
D0 = Pnrm(D0);
                  
% Set up algorithm parameters
lambda = 30e-2;
lambda2 = 12e-2;
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.MaxCGIter = 100;
opt.CGTol = 1e-9;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
opt.sigma = T;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;

close all
lambda2s = [1e-2 5e-2 1e-1 5e-1 1];
%% Dictionary learning
for i = 3%2:numel(sigmas)
    figDir = [topDir,'_sig_',num2str(i)];
    mkdir(figDir)
    % Data
    [y,y_true,N,M,T] = gaus_linear_osc_signal(sigmas(i));
%     plotDataSeq(y_true,topDir,'y_true.gif')
%     for j = 18 %1:numel(lambdas)
    for j = 1:5
        % Solve
        lambda = 30e-2;
        lambda2 = lambda2s(j);
        opt.rho = 50*lambda + 0.5;
        opt.rho2 = 5*lambda2;
        [D,X,Dmin,Xmin,Uvel,Vvel,optinf,obj,relErr] = cbpdndl_cg_OF_multiScales(D0, y, lambda,lambda2s(j), opt, scales);
        
        % Save outputs
        outputs = struct();
        outputs.y = y;
        outputs.D = D;
        outputs.X = X;
        outputs.Dmin = Dmin;
        outputs.Xmin = Xmin;
        outputs.scales = scales;
        outputs.N = N;
        outputs.M = M;
        outputs.T = T;
        outputs.K = K;
        outputs.opt = opt;
        outputs.lambda = lambda;
        outputs.lambda2 = lambda2;
        suffix = sprintf('_sig_%i_lam1_%0.2e_lam2_%0.2e',...
                          j,outputs.lambda,outputs.lambda2);
        save(fullfile(figDir,['output',suffix,'.mat']),'outputs');
        
        % Generate figures
        generateFiguresToy1(figDir,outputs,suffix)

        AD = reSampleCustomArray(N,D,scales);
        ADf = fft2(AD);
        Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));
        plotDataRecon(y,Yhat,figDir,['y_recon',suffix,'.gif'])
        close all
    end
end
%% Load outputs and regen figures with min instead
% for j = 1:3
%     load(fullfile(figDir,sprintf('output_%i.mat',j)))
%     outputs.D = outputs.Dmin;
%     outputs.X = outputs.Xmin;
%     
%     % Generate figures
%     generateFiguresToy1(figDir,outputs,j)
%     
%     AD = reSampleCustomArray(N,Dmin,scales);
%     ADf = fft2(AD);
%     Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(Xmin)),3),'symmetric'));
%     plotDataRecon(y,Yhat,figDir,sprintf('y_recon_%i.gif',j))
%     close all
% end
