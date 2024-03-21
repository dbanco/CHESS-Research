%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = [  5e-2 ];
lambdaHSVals = [1e-8 1e-6 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 0.1 1];
lambdaOFVals = [0    1e-3 2e-3 5e-3 1e-2,...
                2e-2 5e-2 0.1  0.2  0.5,...
                1    2    4    6    8,...
                10   15   20   25   30,...
                17.5 35   40   50   75,...
                100 200 500 1000 2000,...
                1e5];
for j_hs = 6
topDir = ['C:\Users\dpqb1\Documents\Outputs\gaus_example_matched_9_14_23_X0_D0_V0',num2str(lambdaHSVals(j_hs))];
% topDir = '/cluster/home/dbanco02/Outputs/toy1_exp_OF1vel1_matched';

% Experiment Setup
sigmas = 0:0.01:0.05;

% Data parameters
[y,y_true,N,M,T] = gaus_example_matched_multiscale_dl(0.01);
y = reshape(y,[1,N,T]);

% Model Setup
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);
KJ = K*J;
opt = [];
opt.DictFilterSizes = [1;...
                       M];

% Init solution
opt.Y0 = zeros(1,N+M-1,KJ,T);


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
for i = 2%2:numel(sigmas)
    figDir = [topDir,'_sig_',num2str(i)];
    mkdir(figDir)
    
    % Data  
    [y,~,N,M,T] = gaus_example_matched_multiscale_dl();
    y = reshape(y,[1,N,T]);
    center = (M+1)/2;
    
    % Initialization 
    D0(1,round(M/3):round(2*M/3),1) = 1;
    D0 = Pnrm(D0);
    opt.G0 = D0;


    % Rho and sigma params
% opt.rho = 50*lambda + 0.5;
% opt.sigma = T;
    opt.rho = 1e3;%100;
    opt.sigma = 1e3;%100;


    for j_s = 1
        for j_of = 1
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