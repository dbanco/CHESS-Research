%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];

% Experiment Setup
sigmas = 0:0.01:0.1;

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 1000;
opt.NoOFIters = 0;
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
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);

scriptFileName = 'mcdlof_bash.sh';
funcName = 'sim_mcdlof_wrapper';
jobDir = '/cluster/home/dbanco02/jobs/';
k = 1;

datasets = {'pseudo-voigt_unmatched','unmatched'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenter = {0,1};

% Noise level
sig_ind = 1:6;

% SELECTE PARAMETERS:
selected_lam_s = [0,7,11,17,19,26];

% Regularization parameters
% ind1 = 11:15;
% ind2 = [1,30,31];
% ind3 = [5];

ind1 = 1:106;
% ind2 = 2:50;
% ind3 = 2:6;
ind2 = 1;
ind3 = 1;

k = 1;
for s_dataset = 1
dataset = datasets{s_dataset};
for trials = 1
for s_pen = 2
for s_xinit = 1
for s_dinit = 2
for s_dfix = 1
for s_recenter = 2
    opt.Penalty = penalties{s_pen};
    opt.coefInit = xinits{s_xinit};
    opt.dictInit = dinits{s_dinit};
    opt.Dfixed = dfixes{s_dfix};
    opt.Recenter = recenter{s_recenter};
    
    if (opt.Dfixed == 1) && strcmp(opt.dictInit, 'flat')
        continue
    end
   
    topDir = ['/cluster/home/dbanco02/Outputs_6_10_',dataset,'_',opt.Penalty,...
        '_D',opt.dictInit,num2str(opt.Dfixed),...
        '_X',opt.coefInit,num2str(opt.Xfixed),...
        '_recenter',num2str(opt.Recenter),'/results'];
    
    %%% PARAMETER SELECTION SETUP
    % criterion = 'relaxed discrepancy';
    % useMin = 1;
    % relax_param = 1.05;
    % psDir = ['E:\MCDLOF_processing\Outputs_4_19_',dataset,'_',opt.Penalty,...
    % '_D',opt.dictInit,num2str(opt.Dfixed),...
    % '_X',opt.coefInit,num2str(opt.Xfixed),...
    % '_recenter',num2str(opt.Recenter)];
    %%%
    selected_lam_s_inds = [0,9,12,21,30,37];

    for sig_i = sig_ind
        %%% PARAMETER SELECTION
        % inDir = [psDir,'\results_sig_',num2str(sig_i)];
        % [~,y_true,~,~,~] = sim_switch_multiscale_dl(sigmas(sig_i),dataset);
        % [lambda_all,objective] = param_select_3D(inDir,0,criterion,sigmas(sig_i),y_true,useMin,relax_param);
        % j_s_select = find(lambdaVals == lambda_all(1));
        %%%
        j_s_select = selected_lam_s_inds(sig_i);
        for j_s = ind1%j_s_select-2:j_s_select+2
        % for j_s = ind1
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