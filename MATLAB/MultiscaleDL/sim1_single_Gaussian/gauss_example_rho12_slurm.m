%% Do grid search over rho parameters
lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
lambdaRhos = [10,33];

rho1Vals = logspace(-2,3,20);
rho2Vals = logspace(-2,3,20);

% Experiment Setup
sigmas = 0:0.01:0.1;

% Data parameters
[y,y_true,N,M,T] = gaus_example_multiscale_dl();
y = reshape(y,[1,N,T]);

SNR = zeros(T);
for t = 1:T
    SNR(t) = norm(y_true(:,t))/norm(y(1,:,t)-y_true(:,t));
end

% Model Setup
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],10,32, 1/10);
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
D0(1,:) = 1;
D0 = Pnrm(D0);
opt.G0 = D0;

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 1000;
opt.MaxCGIter = 100;
opt.CGTol = 1e-6;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
% Rho and sigma params
% opt.rho = 50*lambda + 0.5;
% opt.sigma = T;
% opt.rho = 1e-1;%1000;
% opt.sigma = 1e-1;%1000;
opt.AutoRho = 0;
opt.AutoRhoPeriod = 1;%10
opt.AutoSigma = 0;
opt.AutoSigmaPeriod = 1;%10
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.UpdateVelocity = 1;
opt.HSiters = 100;
opt.useGpu = 0;

k = 1;
scriptFileName = 'mcdlof_bash.sh';
funcName = 'mcdlof_wrapper_sim1_rho12';
jobDir = '/cluster/home/dbanco02/jobs/';
for i = 3
    for sub_j = 1:numel(lambdaRhos)
        for j_hs = 1     
            for j_of = 1
                for j_rho1 = 1:numel(rho1Vals)
                    for j_rho2 = 1:numel(rho2Vals)
                        opt.rho = rho1Vals(j_rho1);
                        opt.sigma = rho2Vals(j_rho2);
                        topDir = ['/cluster/home/dbanco02/Outputs/gaus_example_10_1_24_X0_D0_V0',num2str(lambdaHSVals(j_hs))];
                        figDir = [topDir,'_sig_',num2str(i)];
                        mkdir(figDir)
                        close all
                        
                        varin = {lambdaVals,lambdaOFVals,lambdaHSVals,lambdaRhos(sub_j),j_of,j_hs,j_rho1,j_rho2,sigmas,i,opt,K,scales,topDir};
                        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
                        k = k + 1;
                    end
                end
            end
        end
    end
end
slurm_write_bash(k-1,jobDir,scriptFileName,sprintf('1-%i',k-1))