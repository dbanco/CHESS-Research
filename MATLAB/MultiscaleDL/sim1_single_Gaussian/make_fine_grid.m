%% Stage 2: Fine grid generation based on coarse results
dataset = 'dissertation_adjust2';
SNRs = [20,16,12,8,4];
sigmas = zeros(numel(SNRs),1);
for i = 1:numel(SNRs)
    sigmas(i) = SNRtoSigma(SNRs(i),dataset);
end
for sig_i = 1:numel(SNRs)
    sigma = sigmas(sig_i);
    rng(1);
    [y,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigma,dataset);
    coarseDir = ['/cluster/home/dbanco02/Outputs_9_29_filter3_dissertation_adjust2_log_Dflat0_Xzeros0/results_trial_1_sig_',num2str(sig_i)];
    fineDir = ['/cluster/home/dbanco02/Outputs_9_29_fine_filter3_dissertation_adjust2_log_Dflat0_Xzeros0/results_trial_1_sig_',num2str(sig_i)];
    jobDir = '/cluster/home/dbanco02/jobs/';
    scriptFileName = 'mcdlof_fine_bash.sh';
    funcName = 'sim_mcdl_reg_wrapper';
    K_select = 4;           % number of coarse candidates to refine
    tolerance = 0.2;        % discrepancy ±20%
    nFine = 10;             % points per axis in fine grid
    zoomFactor = 2;         % factor around coarse center
    
    if ~exist(fineJobDir,'dir'), mkdir(fineJobDir); end
    
    % --- Load all coarse metrics ---
    metricFiles = dir(fullfile(coarseDir,'output_*.mat'));
    metricsAll = [];
    lambdaValsAll = [];
    lambda2ValsAll = [];
    relErrAll = [];
    logPenAll = [];
    regPenAll = [];
    
    for f = 1:numel(metricFiles)
        data = load(fullfile(coarseDir,metricFiles(f).name));
        m = data.outputs.metrics;
        if f == 1
            opt = data.outputs.opt;
        end
        lambdaValsAll(end+1) = m.lambda;
        lambda2ValsAll(end+1)= m.lambda2;
        relErrAll(end+1) = m.rel_error;
        logPenAll(end+1) = m.log_penalty;
        regPenAll(end+1) = m.reg_penalty;
    end
    
    % --- Apply discrepancy principle ---
    target_rel_err = sigma.^2 .* numel(y) ./ sum(y.^2,'all'); % adjust as needed
    withinBand = abs(relErrAll - target_rel_err) ./ target_rel_err <= tolerance;
    
    if ~any(withinBand)
        % pick closest if no run hits band
        [~,bestIdx] = min(abs(relErrAll - target_rel_err));
        candidateIdx = bestIdx;
    else
        candidateIdx = find(withinBand);
    end
    
    % --- Score candidates: normalized combination of rel_error, log_penalty, reg_penalty ---
    v_err = relErrAll(candidateIdx); v_err = (v_err - min(v_err)) / (max(v_err)-min(v_err)+eps);
    v_log = logPenAll(candidateIdx); v_log = (v_log - min(v_log)) / (max(v_log)-min(v_log)+eps);
    v_reg = regPenAll(candidateIdx); v_reg = (v_reg - min(v_reg)) / (max(v_reg)-min(v_reg)+eps);
    
    w_err = 0.5; w_log = 0.25; w_reg = 0.25;
    score = w_err*v_err + w_log*v_log + w_reg*v_reg;
    
    [~,ord] = sort(score);
    K = min(K_select, numel(ord));
    selectedIdx = candidateIdx(ord(1:K));
    
    fprintf('Selected %d coarse candidates for fine grid\n',K);
    
    % --- Generate fine grids around selected parameters ---
    k = 1;
    for i = 1:K
        idx = selectedIdx(i);
        lam = lambdaValsAll(idx);
        lam2 = lambda2ValsAll(idx);
        
        fineLambda = logspace(log10(lam/zoomFactor), log10(lam*zoomFactor), nFine);
        fineLambda2= logspace(log10(lam2/zoomFactor), log10(lam2*zoomFactor), nFine);
        
        for j_s = 1:numel(fineLambda)
            for j_reg = 1:numel(fineLambda2)
                opt.lambda = fineLambda(j_s);
                opt.lambda2 = fineLambda(j_reg);
                varin = {fineLambda,fineLambda2,j_s,j_reg,sigmas,sig_i,opt,fineDir,dataset};
                save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
                k = k + 1;
            end
        end
    end
end

% --- Write SLURM job script ---
slurm_write_bash(k-1,fineJobDir,scriptFileName,sprintf('1-%i',k-1));
fprintf('Wrote %d fine jobs to %s\n', k-1, fineJobDir);
