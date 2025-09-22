%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = logspace(-2,0,150);
lambdaRegVals = [0 logspace(-2,2,99)];

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 300;
opt.MaxCGIter = 100;
opt.CGTol = 1e-6;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
% Rho1 and Rho2 params  
opt.rho1 = 100;
opt.rho2 = 5;
opt.AutoRho1 = 1;
opt.AutoRho1Period = 1;
opt.AutoRho2 = 1;
opt.AutoRho2Period = 1;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.a_min = 1e-4;
opt.lambda_min = 1e-4;
opt.adapt_a = false;
opt.adapt_lambda = false;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.useGpu = false;
opt.Xfixed = 0;
opt.Dfixed = 0;
opt.Recenter = 0;
opt.a = 1;
opt.useMin = false;
opt.AdaptIters = 100;
opt.a_via_lam = true;
opt.l1_iters = 10;
opt.mcdl_init = false;
opt.ism_init = true;
opt.L = 1;

% Multiscale dictionary setup
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);

scriptFileName = 'mcdlof_bash.sh';
funcName = 'sim_mcdl_reg_wrapper';
jobDir = '/cluster/home/dbanco02/jobs/';
k = 1;

datasets = {'pseudo-voigt_unmatched'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenter = {0,1};

% Noise level
dataset = datasets{1};
sig_ind = 1:5;
SNRs = [20,16,12,8,4];
sigmas = zeros(numel(SNRs),1);
for i = sig_ind
    sigmas(i) = SNRtoSigma(SNRs(i),dataset);
end

ind1 = 1:numel(lambdaVals);

% Regularizer
opt.regularizer = 'filter3';

% --- Dataset, Initialization, Parameters ---
for trial = 1
for s_pen = 2
for s_xinit = 1
for s_dinit = 3
for s_dfix = 2
for s_recenter = 1
    opt.Penalty = penalties{s_pen};
    opt.coefInit = xinits{s_xinit};
    opt.dictInit = dinits{s_dinit};
    opt.Dfixed = dfixes{s_dfix};
    opt.Recenter = recenter{s_recenter};
    
    if (opt.Dfixed == 1) && strcmp(opt.dictInit, 'flat')
        continue
    end
   
    topDir = ['/cluster/home/dbanco02/Outputs_9_22_',opt.regularizer,'_',dataset,'_',opt.Penalty,...
        '_D',opt.dictInit,num2str(opt.Dfixed),...
        '_X',opt.coefInit,num2str(opt.Xfixed),...
        '_recenter',num2str(opt.Recenter),'/results_trial_',num2str(trial)];

    for sig_i = 3:5
        for j_s = 90
            for j_reg = [1,5:5:55]
                sim_mcdl_reg_wrapper(lambdaVals,lambdaRegVals,...
                        j_s,j_reg,sigmas,sig_i,opt,topDir,dataset,K,scales);
            end
        end
    end
end
end
end
end
end
end

slurm_write_bash(k-1,jobDir,scriptFileName,sprintf('1-%i',k-1))