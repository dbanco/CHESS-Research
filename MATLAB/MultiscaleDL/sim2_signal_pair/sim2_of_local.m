%% Multiscale 1D dictionary learning toy problem
% Directory
% lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaVals = logspace(-3,1,120);
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,19),10,100];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3 5e-3];

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
opt.rho = 1000;%300;
opt.sigma = 1;%50;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.SigmaScaling = 1.1;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;
opt.SigmaScaling = 1.1;
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
opt.Dfixed = 1;
opt.Recenter = 1;
% opt.Penalty = 'l1-norm';
% opt.Penalty = 'log';
% opt.coefInit = 'zeros';
% opt.dictInit = 'flat';
opt.a = 1;
opt.useMin = false;
opt.AdaptIters = 100;

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

datasets = {'gaussian_tooth_matched_long2',...
            'voigt_tooth_matched_long3',...
            'dissertation_adjust2',...
            'voigt_example_6peaks'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true','mcdl'};
dfixes = {0,1};
recenters = {0,1};

sig_ind = 1:6;

% ind1 = 64%1:numel(lambdaVals);
% ind2 = 15%2:20;
% ind3 = 2%2:6;
% ind2 = 1;
% ind3 = 1;

% --- Dataset, Initialization, Parameters ---
for s0 = 3
dataset = datasets{s0};
for trial = 1
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
    
    topDir = ['E:\MCDLOF_processing\\Outputs_8_19_local_',dataset,'_',opt.Penalty,...
        '_D',opt.dictInit,num2str(opt.Dfixed),...
        '_X',opt.coefInit,num2str(opt.Xfixed),...
        '_recenter',num2str(opt.Recenter),'/results_trial_',num2str(trial)];
        
    selected_lam_s_inds = [0,41,51,55,61,66];
    % selected_lam_s_inds = [43,43,56,58,66,70];
    % selected_lam_of_inds = [3,3,5,13,7,20];
    % selected_lam_hs_inds = [6,5,3,4,2,6];

        % --- Noise level, regularization parameters ---
        for sig_i = 2
        j_s_select = selected_lam_s_inds(sig_i);
        initDir = ['E:\MCDLOF_processing\Outputs_7_23a1_',dataset,'_',opt.Penalty,...
        '_Dflat0','_Xzeros0','_recenter0','\results_trial_1_sig_',num2str(sig_i)];
        % j_of_select = selected_lam_of_inds(sig_i);
        % j_hs_select = selected_lam_hs_inds(sig_i);
        files = dir(fullfile(initDir,['output_j',num2str(j_s_select),'_1_1*.mat']));
        opt.mcdl_file = fullfile(initDir,files(1).name);
        for j_s = j_s_select
        for j_of = 10
        for j_hs = 5
            sim_mcdlof_wrapper3(lambdaVals,lambdaOFVals,lambdaHSVals,...
                j_s,j_of,j_hs,sigmas,sig_i,opt,topDir,dataset,K,scales);
        end
        end
        end
        end
end
end
end
end
end
end
end