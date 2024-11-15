%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];

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

% Initialization 
D0(1,:) = 1;
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
opt.rho = 1e-1;
opt.sigma = 1e-1;
% opt.rho = 50*lambda + 0.5;
% opt.sigma = T;
% opt.rho = 1e3;%100;
% opt.sigma = 1e3;%100;
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
opt.useGpu = 0;

k = 1;
scriptFileName = 'mcdlof_bash.sh';
funcName = 'mcdlof_wrapper_sim1';
jobDir = '/cluster/home/dbanco02/jobs/';
tradeoff = 0.5;
scaleP = [0.4,5.4,9.48,116];
for i = 2:numel(sigmas)
    for j_hs = 1:4     
        topDir = ['/cluster/home/dbanco02/Outputs/gaus_example_9_12_24_coupled',num2str(lambdaHSVals(j_hs))];
        inDir = ['/cluster/home/dbanco02/Outputs/gaus_example_8_28_24_X0_D0_V00_sig_',num2str(i)];
        [lambda_s_sel,j_s] = param_select_lambda_s(inDir,tradeoff,scaleP);
        for j_of = 1:numel(lambdaOFVals)             
            varin = {lambdaVals,lambdaOFVals,lambdaHSVals,j_s,j_of,j_hs,sigmas,i,opt,K,scales,topDir};
            save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
            k = k + 1;        
        end
    end
end
slurm_write_bash(k-1,jobDir,scriptFileName,sprintf('1-%i',k-1))