%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];

% Experiment Setup
sigmas = 0:0.01:0.1;

% Data parameters
[y,y_true,N,M,T] = gaus_example_matched_multiscale_dl();
y = reshape(y,[1,N,T]);

SNR = zeros(T);
for t = 1:T
    SNR(t) = norm(y_true(:,t))/norm(y(1,:,t)-y_true(:,t));
end

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
opt.rho = 300;
opt.sigma = 50;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;
opt.XRelaxParam = 1; %1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.UpdateVelocity = 1;
opt.HSiters = 100;
opt.useGpu = 0;
opt.Xfixed = 0;
opt.Dfixed = 1;

k = 1;
scriptFileName = 'mcdlof_bash.sh';
jobDir = '/cluster/home/dbanco02/jobs/';
for i = 1
    for j_s = 1:5
        for j_hs = 1     
            for j_of = 1
                topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_11_Dfixed_Xzero_A\matched_results';
                figDir = [topDir,'_sig_',num2str(i)];
                mkdir(figDir)
                
                mcdlof_wrapper_sim1_matched_true(lambdaVals,lambdaOFVals,lambdaHSVals,...
                                            j_s,j_of,j_hs,sigmas,i,opt,K,scales,topDir);
            end
        end
    end
end
