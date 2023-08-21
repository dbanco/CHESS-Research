%% Multiscale 1D dictionary learning toy problem
% Directory
% topDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF_Jterm_DXVtrue_150_1';
topDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF_Jterm_Dtrue_X0_Vtrue_150_6noJ_alpha1';
% topDir = '/cluster/home/dbanco02/Outputs/toy1_exp_OF1vel1_matched';
% mkdir(topDir)

% Experiment Setup
sigmas = 0:0.01:0.05;

% Data parameters
[y,~,N,M,T] = gaus_linear_osc_signal_matched(0);
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
% opt.U0 = zeros(1,N,KM,T);

% Init dictionary
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = zeros(1,M,K); 

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.MaxCGIter = 100;
opt.CGTol = 1e-8;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-8;
opt.sigma = T;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.UpdateVelocity = 0;
opt.Smoothness = 1;%opt.Smoothness = 1e-8;
opt. HSiters = 300;

close all
lambdas = [1e-2 10e-2];
lambda2s = [5 10 20 100];
%% Dictionary learning
for i = 1%2:numel(sigmas)
    figDir = [topDir,'_sig_',num2str(i)];
    mkdir(figDir)
    % Data
    [y,y_true,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched(sigmas(i));
    [u,v,Fy,Fx,Ft]  = computeHornSchunkDict(Xtrue,K,opt.Smoothness,opt.HSiters);
    Jterm = Ft == 0;
%     D0(1,:,1) = Dtrue(:,:,1);
%     D0(1,:,2) = Dtrue(:,:,2);
    D0(1,round(M/3):round(2*M/3),1) = 1;
    D0(1,round(M/4):round(3*M/4),2) = 1;
    D0 = Pnrm(D0);
    Jof = sum(vec(opticalFlowOp(Xtrue,u,v,K,0,Jterm)).^2);
%     plotOpticalFlow(Xtrue,K)
%     plotDataSeq(y_true,topDir,'y_true.gif')
%     for j = 18 %1:numel(lambdas)
    for j = 1:4
        % Solve
%         lambda = 1e-2;
%         lambda2 = 1;
        lambda = 20e-2;
        lambda2 = lambda2s(j);
        opt.rho = 50*lambda + 0.5;
        opt.rho2 = 50*lambda2;
%         opt.Y0 = Xtrue;
        opt.Y0 = zeros(1,N,KM,T);
%         opt.Y0(:,:,1,:) = y;
%         opt.Y0(:,:,2,:) = y;
%         tic
        Jterm = [];
        [D1,X1,Dindep,Xindep,~,~,~,~,~] = cbpdndl_cg_OF_multiScales_gpu(Dtrue, y, lambda,0, opt, scales,u,v,Jterm);
        opt.Y0 = Xindep;
        [D,X,Dmin,Xmin,Uvel,Vvel,optinf,obj,relErr] = cbpdndl_cg_OF_multiScales_gpu(Dtrue, y, lambda,lambda2, opt, scales,u,v,Jterm);
%         toc
%         tic
%         [D,X,Dmin,Xmin,Uvel,Vvel,optinf,obj,relErr] = cbpdndl_cg_OF_multiScales(Dtrue, y, lambda,lambda2, opt, scales,u,v,Jterm);
%         toc
%         [D, X, optinf, obj, relErr] = cbpdndl_cg_multiScales(D0, y, lambda, opt, scales);
        
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
        suffix = sprintf('_j%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e',...
                          j,sigmas(i),outputs.lambda,outputs.lambda2);
        save(fullfile(figDir,['output',suffix,'.mat']),'outputs');
        
        % Generate figures
        indepOutputs = outputs;
        indepOutputs.D = Dindep;
        indepOutputs.X = Xindep;
        indepOutputs.Dmin = Dindep;
        indepOutputs.Xmin = Xindep;
        indepOutputs.lambda2 = 0;
        suffixIndep = sprintf('_indep_j%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e',...
                          j,sigmas(i),outputs.lambda,outputs.lambda2);
        generateFiguresToy1(figDir,indepOutputs,suffixIndep)

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
