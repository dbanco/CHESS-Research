%% Multiscale 1D dictionary learning toy problem
% Directory
% lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaVals = logspace(-3,1,120);
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
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
opt.rho = 300;
opt.sigma = 50;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.UpdateVelocity = 1;
opt.HSiters = 100;
opt.useGpu = 0;
opt.Xfixed = 0;
opt.Dfixed = 0;
opt.Recenter = 0;
opt.a = 1;

% Multiscale dictionary setup
scriptFileName = 'mcdlof_bash.sh';
funcName = 'sim_mcdlof_wrapper3_mmpad';
jobDir = '/cluster/home/dbanco02/jobs/';
k = 1;

datasets = {'mmpad_ring1'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenters = {0,1};

sig_ind = 1;
ind1 = 1:numel(lambdaVals);
ind2 = 1;%2:50;
ind3 = 1;%2:6;

% SELECTED PARAMS ind1: 0,4,4,9
% selected_lam_s = [0,6,6,9];
% j_s_select = find(lambdaVals == selected_lam_s_vec(sig_i));

% --- Dataset, Initialization, Parameters ---
for s0 = 1
dataset = datasets{s0};
for trials = 1
for s_pen = 2
for s_xinit = 1
for s_dinit = 2
for s_dfix = 1
for s_recenter = 1
for K = 1:6
    opt.Penalty = penalties{s_pen};
    opt.coefInit = xinits{s_xinit};
    opt.dictInit = dinits{s_dinit};
    opt.Dfixed = dfixes{s_dfix};
    opt.Recenter = recenters{s_recenter};

    scales = cell(K,1);
    for i = 1:K
        scales{i} = genRationals([0;1],[1;1],8,8, 1/6);
    end
    J = size(scales{1},2);

    if (opt.Dfixed == 1) && strcmp(opt.dictInit, 'flat')
        continue
    end

    topDir = ['/cluster/home/dbanco02/Outputs_7_2_',dataset,'_',opt.Penalty,...
        '_D',opt.dictInit,num2str(opt.Dfixed),...
        '_X',opt.coefInit,num2str(opt.Xfixed),...
        '_recenter',num2str(opt.Recenter),'/results'];
        
        % --- Regularization parameters ---
            % ind1 = selected_lam_s(sig_i);
        for j_s = ind1
        for j_of = ind2
        for j_hs = ind3
            
            varin = {lambdaVals,lambdaOFVals,lambdaHSVals,...
                    j_s,j_of,j_hs,opt,topDir,dataset,K,scales};
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