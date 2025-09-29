%% Stage 2: Fine grid search around best coarse runs
clear; clc;

% --- Directories ---
coarseDir = '/cluster/home/dbanco02/Outputs_9_29_filter3_dissertation_adjust2_log_Dflat0_Xzeros0/results_trial_1';
jobDir    = '/cluster/home/dbanco02/jobs_fine/';
if ~exist(jobDir,'dir'); mkdir(jobDir); end

scriptFileName = 'mcdlof_bash_fine.sh';
funcName = 'sim_mcdl_reg_wrapper';

% --- Load coarse results ---
files = dir(fullfile(coarseDir,'metrics_*.mat'));

allRelErrors = [];
paramList = [];

for f = 1:numel(files)
    data = load(fullfile(coarseDir,files(f).name));
    if ~isfield(data,'metrics'); continue; end
    m = data.metrics;

    allRelErrors(end+1) = m.error;
    paramList(end+1,:)  = [m.lambda, m.lambda2];
end

% --- Pick best candidates ---
[~,bestIdx] = mink(allRelErrors,5);   % top 5 by relative error
bestParams = paramList(bestIdx,:);

% --- Create fine grids around best params ---
fineLambdaVals = [];
fineLambdaRegVals = [];
for i = 1:size(bestParams,1)
    lam    = bestParams(i,1);
    lamReg = bestParams(i,2);

    % zoom ±2× around the coarse point, 10 log-spaced values
    fineLambdaVals     = [fineLambdaVals     logspace(log10(lam/2), log10(lam*2), 10)];
    fineLambdaRegVals  = [fineLambdaRegVals  logspace(log10(max(lamReg/2,1e-4)), log10(max(lamReg*2,1e-4)), 10)];
end

fineLambdaVals    = unique(fineLambdaVals);
fineLambdaRegVals = unique(fineLambdaRegVals);

% --- Set up fine jobs ---
k = 1;
opt = getDefaultOptions(); % your helper with standard options

dataset = 'dissertation_adjust2';
K = 2;
scales = {genRationals([0;1],[1;1],8,8,1/6), genRationals([0;1],[1;1],8,8,1/6)};
sigmas = ... % same as stage 1

for lam = fineLambdaVals
    for lamReg = fineLambdaRegVals
        varin = {lam, lamReg, sigmas, opt, dataset, K, scales};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

slurm_write_bash(k-1,jobDir,scriptFileName,sprintf('1-%i',k-1));
