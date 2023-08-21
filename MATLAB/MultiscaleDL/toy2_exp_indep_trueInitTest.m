%% Multiscale 1D dictionary learning toy problem
% Directory
% topDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF_Jterm_DXVtrue_150_1';
topDir = 'C:\Users\dpqb1\Documents\Outputs\toy2_exp_optFlow7_20_Xzero_D0_indep1';
% topDir = '/cluster/home/dbanco02/Outputs/toy1_exp_OF1vel1_matched';

% Experiment Setup
sigmas = 0:0.01:0.05;

% Data parameters
[y,~,K,J,N,M,T,~,~,scales] = gaus_linear_osc_signal_matched_small(0);
y = reshape(y,[1,N,T]);

% Model Setup
KJ = K*J;
opt = [];
opt.DictFilterSizes = [1,1;...
                       M,M];

% Init solution
opt.Y0 = zeros(1,N,KJ,T);
% opt.U0 = zeros(1,N,KJ,T);

% Init dictionary
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = zeros(1,M,K); 

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 5000;
opt.MaxCGIter = 100;
opt.CGTol = 1e-6;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;%10
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;%10
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.UpdateVelocity = 1;
opt.Smoothness = 10;%1e-6;%opt.Smoothness = 1e-8;
opt. HSiters = 100;

close all
lambdas = [1e-2 10e-2];
lambda2s = [1e-3 1e-2 0.1 1 5 10 20 100 1000];
%% Dictionary learning
for i = 2%2:numel(sigmas)
    figDir = [topDir,'_sig_',num2str(i)];
    mkdir(figDir)
    % Data  
    [y,y_true,K,J,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small(sigmas(i));
    norm(y(:)-y_true(:))/norm(y_true(:))
    % Independent solution initialization 
    D0(1,round(M/3):round(2*M/3),1) = 1;
    D0(1,round(M/4):round(3*M/4),2) = 1;
    D0 = Pnrm(D0);
    opt.Y0 = Xtrue;
%     opt.Y0 = zeros(size(Xtrue));
    opt.G0 = Dtrue;
    lambda = 1e-3;
    lambda2 = 0;
% opt.rho = 50*lambda + 0.5;
% opt.sigma = T;
    opt.rho = 1e2;
    opt.sigma = 1e2;
    opt.UpdateVelocity = 0;
%     load('indepOutputs_7_7.mat')

    [u,v,~,~,~]  = computeHornSchunkDictPaperLS(Xtrue,K,[],[],opt.Smoothness,opt.HSiters);
    [Dindep,Xindep,D1,X1,~,~,optinf,~,~] = cbpdndl_cg_OF_multiScales_gpu(Dtrue, y, lambda,lambda2, opt, scales,u,v);
    [uIndep,vIndep,~,~,~]  = computeHornSchunkDictPaperLS(Xindep,K,[],[],opt.Smoothness,opt.HSiters);
    indepOutputs.y = y;
    indepOutputs.D = Dindep;
    indepOutputs.X = Xindep;
    indepOutputs.Dmin = D1;
    indepOutputs.Xmin = X1;
    indepOutputs.scales = scales;
    indepOutputs.N = N;
    indepOutputs.M = M;
    indepOutputs.T = T;
    indepOutputs.K = K;
    indepOutputs.opt = opt;
    indepOutputs.lambda = lambda;
    indepOutputs.lambda2 = lambda2;
    indepOutputs.Uvel = uIndep;
    indepOutputs.Vvel = vIndep;
    suffixIndep = sprintf('_indep_sig_%0.2e_lam1_%0.2e',...
                      sigmas(i),indepOutputs.lambda);
    generateFiguresToy1(figDir,indepOutputs,suffixIndep,[4,4])
    opt.Y0 = indepOutputs.X;
    opt.Y0 = Xindep;
    opt.rho = optinf.rho;
    opt.sigma = optinf.sigma;

%     plotOpticalFlow2(Xtrue,K,opt)
%     plotDataSeq(y_true,topDir,'y_true.gif')
%     for j = 18 %1:numel(lambdas)

%     for j = 1:3
%         % Optical flow coupled solution
%         lambda2 = lambda2s(j);     
%         opt.UpdateVelocity = 1;
%         [D,X,Dmin,Xmin,Uvel,Vvel,optinf,obj,relErr] = cbpdndl_cg_OF_multiScales_gpu(Dindep, y, lambda,lambda2, opt, scales,u,v);
%         
%         % Save outputs
%         outputs = struct();
%         outputs.y = y;
%         outputs.D = D;
%         outputs.X = X;
%         outputs.Dmin = Dmin;
%         outputs.Xmin = Xmin;
%         outputs.scales = scales;
%         outputs.N = N;
%         outputs.M = M;
%         outputs.T = T;
%         outputs.K = K;
%         outputs.opt = opt;
%         outputs.lambda = lambda;
%         outputs.lambda2 = lambda2;
%         outputs.Uvel = Uvel;
%         outputs.Vvel = Vvel;
%         suffix = sprintf('_j%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e',...
%                           j,sigmas(i),outputs.lambda,outputs.lambda2);
%         save(fullfile(figDir,['output',suffix,'.mat']),'outputs');
%         
%         % Generate figures
%         generateFiguresToy1(figDir,outputs,suffix)
% %         generateFiguresToy1min([figDir,'min'],outputs,suffix)
% %         generateFiguresToy1([figDir,'indep'],inde,suffix)
% 
%         AD = reSampleCustomArray(N,D,scales);
%         ADf = fft2(AD);
%         Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));
%         plotDataRecon(y,Yhat,figDir,['y_recon',suffix,'.gif'])
%         close all
%     end
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

%% View optical flow
% window =1:300;
% X = outputs.X;
% Uvel = outputs.Uvel;
% Vvel = outputs.Vvel;
% 
% T = 10;
% window =140:200;
% J = 22;
% fig = figure;
% fig.Position = [1 1 1.6452e+03 554.8000];
% for t = 1:T
%     framePlot = squeeze(X(1,window,1:J,t));
%     subplot(3,1,1)
%     imagesc(framePlot')
%     hold on
%     quiver(Vvel(window,1:J,t)',Uvel(window,1:J,t)')
%     q = findobj(gca,'type','Quiver');
%     q.Color = 'w';
%     hold off

%     framePlot = squeeze(Xindep(1,window,1:J,t));
%     subplot(3,1,2)
%     imagesc(framePlot')
%     hold on
%     quiver(vIndep(window,1:J,t)',uIndep(window,1:J,t)')
%     q = findobj(gca,'type','Quiver');
%     q.Color = 'w';
%     hold off
% 
%     framePlot = squeeze(Xtrue(1,window,1:J,t));
%     subplot(3,1,3)
%     imagesc(framePlot')
%     hold on
%     quiver(v(window,1:J,t)',u(window,1:J,t)')
%     q = findobj(gca,'type','Quiver');
%     q.Color = 'w';
%     hold off

%     pause()
% end

