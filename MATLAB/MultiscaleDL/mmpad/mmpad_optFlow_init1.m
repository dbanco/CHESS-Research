%% Multiscale 1D dictionary learning MMPAD

% Load MMPAD data
y = loadMMPAD1D(1,1);
[N,T] = size(y);
y = reshape(y,[1,N,T]);

% Model Setup
K = 2;
M = 75;
center = (M+1)/2;
denLim = 8;
scales{1} = genRationals([0;1],[1;1],denLim,denLim, 1/6);
scales{2} = genRationals([0;1],[1;1],denLim,denLim, 1/6);
J = size(scales{1},2);
KJ = K*J;
opt = [];
opt.DictFilterSizes = [1,1;...
                       M,M];

Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

% Init solution

% load("C:\Users\dpqb1\Documents\Outputs\mmpad_ring1_optFlow\output_j_5_10_5_lam1_2.00e-02_lam2_5.00e-01_lam3_1.00e-03.mat")
% D0 = outputs.D;
% Uvel = outputs.Uvel;
% Vvel = outputs.Vvel;
% opt.Y0 = outputs.X;
opt.Y0 = zeros(1,N+M-1,KJ,T);
Uvel = squeeze(zeros(size(opt.Y0)));
Vvel = squeeze(zeros(size(opt.Y0)));

D0 = zeros(1,M,K); 
D0(1,round(M/3):round(2*M/3),1) = 1;
D0(1,round(M/4):round(3*M/4),2) = 1;
D0 = Pnrm(D0);

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
opt.HSiters = 100;
opt.rho = 1e2;%100;
opt.sigma = 1e2;%100;
opt.G0 = D0;

% Test parameters
lambdaVals = [  1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1 1 2 5];
lambdaHSVals = [1e-8 1e-6 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 0.1 1];
lambdaOFVals = [0    1e-3 2e-3 5e-3 1e-2,...
                2e-2 5e-2 0.1  0.2  0.5,...
                1    2    4    6    8,...
                10   15   20   25   30,...
                17.5 35   40   50   75,...
                100 200 500 1000 2000,...
                1e5];

% Output Dir
outDir = 'C:\Users\dpqb1\Documents\Outputs\mmpad_ring1_optFlow_init0';
mkdir(outDir)


for j_hs = 5
    for j_s = 5
        for j_of = [8,12,16,24]
            % Optical flow coupled solution
            lambda = lambdaVals(j_s);
            lambda2 = lambdaOFVals(j_of);
            lambda3 = lambdaHSVals(j_hs);
            opt.Smoothness = lambda3;
            
            [D,X,Dmin,Xmin,Uvel,Vvel,optinf,obj,relErr] = cbpdndl_cg_OF_multiScales_gpu_zpad_center(D0, y, lambda,lambda2, opt, scales,Uvel,Vvel);
            
            % Save outputs
            outputs = saveOutputs(y,D,X,Dmin,Xmin,Uvel,Vvel,lambda,lambda2,lambda3,N,M,K,T,scales,center,opt,j_s,j_of,j_hs,outDir);
        end
    end
end
