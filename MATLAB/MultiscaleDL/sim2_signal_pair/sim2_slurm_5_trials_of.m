%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = logspace(-2,0,150);
lambdaOFVals = [0 logspace(-4,0,20)];
lambdaHSVals = [0 logspace(-4,0,10)];

% Experiment Setup
sigmas = 0:0.01:0.1;

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 1000;
opt.MaxCGIter = 100;
opt.NoOFIters = 0;
opt.CGTol = 1e-6;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
% Rho and sigma params
opt.rho = 1000;%300
opt.sigma = 500;%50
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.a_min = 1e-4;
opt.lambda_min = 1e-4;
opt.adapt_a = true;
opt.adapt_lambda = false;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.UpdateVelocity = 1;
opt.HSiters = 100;
opt.useGpu = 0;
opt.Xfixed = 0;
opt.Dfixed = 0;
opt.Recenter = 0;
opt.a = 1;
opt.useMin = false;
opt.AdaptIters = 100;
opt.a_via_lam = true;
opt.l1_iters = 10;

% Multiscale dictionary setup
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],8,8, 1/6);
scales{2} = genRationals([0;1],[1;1],8,8, 1/6);
J = size(scales{1},2);

scriptFileName = 'mcdlof_bash.sh';
funcName = 'sim_mcdlof_wrapper3';
jobDir = '/cluster/home/dbanco02/jobs/';
k = 1;

datasets = {'dissertation_adjust2'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true','mcdl'};
dfixes = {0,1};
recenters = {0,1};

sig_ind = 1:6;

selected_lam_s_inds = [9,39,65,78,95,96];
selected_lam_of_inds = [10,21,14,16,12,12];
selected_lam_hs_inds = [5,3,3,6,9,6];


% --- Dataset, Initialization, Parameters ---
for s0 = 1
dataset = datasets{s0};
for trial = 1:50
for s_pen = 2
for s_xinit = 1
for s_dinit = 4
for s_dfix = 1
for s_recenter = 1
    opt.Penalty = penalties{s_pen};
    opt.coefInit = xinits{s_xinit};
    opt.dictInit = dinits{s_dinit};
    opt.Dfixed = dfixes{s_dfix};
    opt.Recenter = recenters{s_recenter};

    if (opt.Dfixed == 1) && strcmp(opt.dictInit, 'flat')
        continue
    end
    
    topDir = ['/cluster/home/dbanco02/Outputs_8_25_of_',dataset,'_',opt.Penalty,...
        '_D',opt.dictInit,num2str(opt.Dfixed),...
        '_X',opt.coefInit,num2str(opt.Xfixed),...
        '_recenter',num2str(opt.Recenter),'/results_trial_',num2str(trial)];

    % --- Noise level, regularization parameters ---
    for sig_i = sig_ind
        j_s_select = selected_lam_s_inds(sig_i);
        j_of_select = selected_lam_of_inds(sig_i);
        j_hs_select = selected_lam_hs_inds(sig_i);

        initDir = ['/cluster/home/dbanco02/Outputs_8_25_indep_trials_',dataset,'_',opt.Penalty,...
                '_Dflat0','_Xzeros0','_recenter0','/results_trial_',num2str(trial),'_sig_',num2str(sig_i)];
    
        files = dir(fullfile(initDir,['output_j',num2str(j_s_select),'_1_1*.mat']));
        opt.mcdl_file = fullfile(initDir,files(1).name);

        j_s = j_s_select;
        j_of = j_of_select;
        j_hs = j_hs_select;
        varin = {lambdaVals,lambdaOFVals,lambdaHSVals,...
                j_s,j_of,j_hs,sigmas,sig_i,opt,topDir,dataset,K,scales};
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
slurm_write_bash(k-1,jobDir,scriptFileName,sprintf('1-%i',k-1))