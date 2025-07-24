%% Multiscale 1D dictionary learning toy problem
% Directory
% lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaVals = logspace(-3,1,120);
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,5)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 1e-1 1];

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

datasets = {'sim2_gaussian_tooth_matched','sim2_gaussian_tooth_unmatched',...
            'sim2_gaussian_tooth_matched2','sim2_gaussian_tooth_unmatched2',...
            'dissertation','dissertation_adjust2','dissertation_long',...
            'dissertation_long_separate','voigt_tooth_matched',...
            'gaussian_tooth_matched','gaussian_tooth_matched_long',...
            'gaussian_tooth_matched_long2'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenters = {0,1};

sig_ind = 2:6;

% ind1 = 50:2:70;
% ind2 = 2:11;
% ind3 = 2:9;

% ind1 = 30:39;
% ind2 = 2:11;
% ind3 = 2:9;

ind1 = 1:120;
ind2 = 1;
ind3 = 1;

% --- Dataset, Initialization, Parameters ---
for s0 = 6
dataset = datasets{s0};
for trial = 1
for s_pen = 1
for s_xinit = 1
for s_dinit = 2
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
    
    topDir = ['/cluster/home/dbanco02/Outputs_7_23_',dataset,'_',opt.Penalty,...
        '_D',opt.dictInit,num2str(opt.Dfixed),...
        '_X',opt.coefInit,num2str(opt.Xfixed),...
        '_recenter',num2str(opt.Recenter),'/results_trial_',num2str(trial)];
        
    % selected_lam_s_inds = [33,36,51,57,62,66];
    % selected_lam_s_inds = [43,43,56,58,66,70];
    % selected_lam_of_inds = [3,3,5,13,7,20];
    % selected_lam_hs_inds = [6,5,3,4,2,6];

        % --- Noise level, regularization parameters ---
        for sig_i = sig_ind
        % j_s_select = selected_lam_s_inds(sig_i);
        % j_of_select = selected_lam_of_inds(sig_i);
        % j_hs_select = selected_lam_hs_inds(sig_i);

        for j_s = ind1
        for j_of = ind2
        for j_hs = ind3
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
end
end
end

slurm_write_bash(k-1,jobDir,scriptFileName,sprintf('1-%i',k-1))