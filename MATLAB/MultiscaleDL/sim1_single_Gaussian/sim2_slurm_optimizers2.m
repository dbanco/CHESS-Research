%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = logspace(-1.5,1,150);
lambdaRegVals = [0 logspace(-2,2,99)];

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 1000;
opt.MaxCGIter = 100;
opt.CGTol = 1e-10;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-10;
% Rho and sigma params
opt.rho1 = 1000;
opt.rho2 = 500;
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
opt.useMin = true;
opt.AdaptIters = 100;
opt.a_via_lam = true;
opt.l1_iters = 20;
opt.mcdl_init = 100;
opt.ism_init = true;
opt.L = 1;
opt.tau = 1e-2;

scriptFileName = 'mcdlof_bash.sh';
funcName = 'sim_mcdl_reg_wrapper';
jobDir = '/cluster/home/dbanco02/jobs/';
k = 1;

datasets = {'dissertation_adjust2'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenter = {0,1};

% Noise level
dataset = datasets{1};
sig_ind = 1:4;
SNRs = [20,17,14,11,8];
sigmas = zeros(numel(SNRs),1);
for i = sig_ind
    sigmas(i) = SNRtoSigma(SNRs(i),dataset);
end

ind1 = 1:numel(lambdaVals);

% Regularizer
opt.regularizer = 'softmin';
optimizers = {'LBFGS', 'LinMM', 'QuadMM'};

% selected_lam = [14,26,39,50,63];

% --- Dataset, Initialization, Parameters ---
for trial = 1
for s_pen = 2
for s_xinit = 1
for s_dinit = 2
for s_dfix = 1
for s_optim = 1
    opt.Penalty = penalties{s_pen};
    opt.coefInit = xinits{s_xinit};
    opt.dictInit = dinits{s_dinit};
    opt.Dfixed = dfixes{s_dfix};
    opt.optimizer = optimizers{s_optim};
    
    if (opt.Dfixed == 1) && strcmp(opt.dictInit, 'flat')
        continue
    end
   
    topDir = ['/cluster/home/dbanco02/Outputs_10_3_reg_',opt.regularizer,'_',opt.optimizer,...
        '_',dataset,'_',opt.Penalty,...
        '_D',opt.dictInit,num2str(opt.Dfixed),...
        '_X',opt.coefInit,num2str(opt.Xfixed),...
        '/results_trial_',num2str(trial)];

    for sig_i = 1:4
        jj = selected_lam(sig_i);
        for j_s = jj:(jj+1)
            for j_reg = 2:100
                varin = {lambdaVals,lambdaRegVals,j_s,j_reg,sigmas,...
                         sig_i,opt,topDir,dataset};
                save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
                k = k + 1;
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