%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = [  5e-2 6e-2 7e-2 8e-2 9e-2];
lambdaHSVals = [1e-8 1e-6 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 0.1 1];
lambdaOFVals = [0    1e-3 2e-3 5e-3 1e-2,...
                2e-2 5e-2 0.1  0.2  0.5,...
                1    2    4    6    8,...
                10   15   20   25   30,...
                17.5 35   40   50   75,...
                100 200 500 1000 2000,...
                1e5];
for j_hs = 4
topDir = ['C:\Users\dpqb1\Documents\Outputs\toy4_center_exp_optFlow8_11_24_X0_D0_V0_zpad_HS2',num2str(lambdaHSVals(j_hs))];

% Experiment Setup
sigmas = [0.01,0.05,0.1,0.15,0.2];

% Data parameters
[y,~,K,J,N,M,T,~,~,scales] = gaus_linear_osc_signal_matched_small_zpad3_center(0);
y = reshape(y,[1,N,T]);

% Model Setup
KJ = K*J;
opt = [];
opt.DictFilterSizes = [1,1;...
                       M,M];

% Init solution
opt.Y0 = zeros(1,N+M-1,KJ,T);
% opt.U0 = zeros(1,N,KJ,T);

% Init dictionary
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
D0 = zeros(1,M,K); 

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 500;
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
opt.Smoothness = lambdaHSVals(j_hs);%1e-6;%opt.Smoothness = 1e-8;
opt.HSiters = 100;

close all
%% Dictionary learning
for i = 1:numel(sigmas)
    figDir = [topDir,'_sig_',num2str(i)];
    mkdir(figDir)
    
    % Data  
    [y,y_true,K,J,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad3_center(sigmas(i));
    center = (M+1)/2;
    
    % Initialization 
    D0(1,round(M/3):round(2*M/3),1) = 1;
    D0(1,round(M/4):round(3*M/4),2) = 1;
    D0 = Pnrm(D0);
    opt.Y0 = Xtrue;
    opt.Y0 = zeros(size(Xtrue));
    opt.G0 = D0;
%     load("C:\Users\dpqb1\Documents\Outputs\toy3_center_exp_optFlow8_17_X0_D0_V0_zpad_HS0.002_sig_2\output_j10_sig_1.00e-02_lam1_6.00e-02_lam2_5.00e-01.mat")
%     D0 = outputs.D;
%     opt.Y0 = outputs.X;
%     opt.G0 = D0;

    % Rho and sigma params
% opt.rho = 50*lambda + 0.5;
% opt.sigma = T;
    opt.rho = 1e3;%100;
    opt.sigma = 1e3;%100;
    opt.UpdateVelocity = 0;

    for j_s = 4:5
        for j_of = 8:12
            % Optical flow coupled solution
            lambda = lambdaVals(j_s);
            lambda2 = lambdaOFVals(j_of);
            [Uvel,Vvel,Fx,Fy,Ft] = computeHornSchunkDictPaperLS(opt.Y0,K,[],[],opt.Smoothness/lambda2,opt.HSiters);
            opt.UpdateVelocity = 1;
            [D,X,Dmin,Xmin,Uvel,Vvel,optinf,obj,relErr] = cbpdndl_cg_OF_multiScales_gpu_zpad_center(D0, y, lambda,lambda2, opt, scales,Uvel,Vvel);
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
            outputs.Uvel = Uvel;
            outputs.Vvel = Vvel;
            suffix = sprintf('_j%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e',...
                              j_of,sigmas(i),outputs.lambda,outputs.lambda2);
            save(fullfile(figDir,['output',suffix,'.mat']),'outputs');
            
            % Generate figures
            generateFiguresToy1zpad_center(figDir,outputs,suffix,[4,4]);
    %         generateFiguresToy1min([figDir,'min'],outputs,suffix)
    %         generateFiguresToy1([figDir,'indep'],inde,suffix)
    
            AD = reSampleCustomArrayCenter(N,D,scales,center);
            AD = padarray(AD,[0 M-1 0],0,'post');
            ADf = fft2(AD);
            Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
            plotDataRecon(y,Yhat,figDir,['y_recon',suffix,'.gif'])
            close all
        end
    end
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

%% View optical flow
% 
% X = outputs.X;
% Uvel = outputs.Uvel;
% Vvel = outputs.Vvel;
% [u,v,~,~,~]    = computeHornSchunkDictPaperLS(Xtrue, K,[],[],1,opt.HSiters);
% [u2,v2,~,~,~]  = computeHornSchunkDictPaperLS2(Xtrue,K,[],[],1,opt.HSiters);
% [u3,v3,~,~,~]  = computeHornSchunkDictPaperLS3(Xtrue,K,[],[],1,opt.HSiters);
% T = 20;
% window = 1:55;
% J = 16;
% fig = figure;
% fig.Position = [1 1 1.6452e+03 554.8000];
% 
% 
% for t = 1:T
%     framePlot = squeeze(Xtrue(1,window,1:J,t));
%     subplot(3,1,1)
%     imagesc(framePlot')
%     hold on
%     quiver(v(window,1:J,t)',u(window,1:J,t)')
%     q = findobj(gca,'type','Quiver');
%     q.Color = 'w';
%     hold off
% 
%     framePlot = squeeze(Xtrue(1,window,1:J,t));
%     subplot(3,1,2)
%     imagesc(framePlot')
%     hold on
%     quiver(v2(window,1:J,t)',u2(window,1:J,t)')
%     q = findobj(gca,'type','Quiver');
%     q.Color = 'w';
%     hold off
% 
%     framePlot = squeeze(Xtrue(1,window,1:J,t));
%     subplot(3,1,3)
%     imagesc(framePlot')
%     hold on
%     quiver(v3(window,1:J,t)',u3(window,1:J,t)')
%     q = findobj(gca,'type','Quiver');
%     q.Color = 'w';
%     hold off
% 
%     pause()
% end

