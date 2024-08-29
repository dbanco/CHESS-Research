%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];

% Experiment Setup
sigmas = 0:0.01:0.1;

% Data parameters
[yn,y,N,M,T] = generate_signal_pair();
y = reshape(y,[1,N,T]);

% Model Setup
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],8,8, 1/6);
scales{2} = genRationals([0;1],[1;1],8,8, 1/6);
J = size(scales{1},2);
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
D0(1,round(M/3):round(2*M/3),1) = 1;
D0(1,round(M/4):round(3*M/4),2) = 1;
D0 = Pnrm(D0);
opt.G0 = D0;

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 500;
opt.MaxCGIter = 100;
opt.CGTol = 1e-6;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
% Rho and sigma params
% opt.rho = 50*lambda + 0.5;
% opt.sigma = T;
opt.rho = 1e-1;%1000;
opt.sigma = 1e-1;%1000;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;%10
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;%10
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.UpdateVelocity = 0;
opt.HSiters = 100;
opt.useGpu = 0;

k = 1;
scriptFileName = 'mcdlof_bash.sh';
funcName = 'mcdlof_wrapper_sim2';
jobDir = '/cluster/home/dbanco02/jobs/';

tradeoff = 1.5;
scaleP = [0.157,5.12,19.02,52.406];

close all
%% Dictionary learning
for i = 2:numel(sigmas)
    inDir = ['/cluster/home/dbanco02/Outputs/signal_pair_8_29_24_X0_D0_V00_sig_',num2str(i)];
    [lambda_s_sel,j_s] = param_select_lambda_s(inDir,tradeoff,scaleP);
    for j_hs = 1
        for j_of = 1
            topDir = ['/cluster/home/dbanco02/Outputs/signal_pair_8_29_24_coupled',num2str(lambdaHSVals(j_hs))];
            varin = {lambdaVals,lambdaOFVals,lambdaHSVals,j_s,j_of,j_hs,sigmas,i,opt,K,scales,topDir};
            save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
            k = k + 1; 
        end
    end
end
