% Evaluate Sim 1 Indep and OF
lambdaVals = logspace(-2,0,150);
lambdaOFVals = [0 logspace(-4,0,20)];
lambdaHSVals = [0 logspace(-4,0,10)];

% Experiment Setup
sigmas = 0:0.01:0.1;
sig_ind = 1:6;

selected_lam_s = [90,96,108,112,117,119];
selected_lam_of = [13,10,3,9,10,12];
selected_lam_hs = [2,3,4,11,9,11];

dataset = "pseudo-voigt_unmatched";
topDir1 = "E:\MCDLOF_processing\Outputs_8_19_indep_pseudo-voigt_unmatched_log_Dflat0_Xzeros0_recenter0";
topDir2 = "E:\MCDLOF_processing\Outputs_8_19_of_pseudo-voigt_unmatched_log_Dmcdl0_Xzeros0_recenter0";

outTable = zeros(8,6,2);

% Indep Metrics
for sig_i = sig_ind
    j_s = selected_lam_s(sig_i);
    j_of = selected_lam_of(sig_i);
    j_hs = selected_lam_hs(sig_i);

    lambda1 = lambdaVals(j_s);
    lambda2 = lambdaOFVals(j_of);
    lambda3 = lambdaHSVals(j_hs);
    HSiters = 100;

    outFile = fullfile(topDir1,...
        ['results_trial_1_sig_',num2str(sig_i)],...
        ['output_j',num2str(j_s),'_1_1*.mat']);
    dFiles = dir(outFile);
    loadFile = fullfile(dFiles(1).folder,dFiles(1).name);
    load(loadFile)

    [outVector, outLabels] = compute_recov_metrics(outputs,dataset,sigmas(sig_i),...
        lambda1,lambda2,lambda3,HSiters);
    outTable(:,sig_i,1) = outVector;
end

% OF Metrics
for sig_i = sig_ind
    j_s = selected_lam_s(sig_i);
    j_of = selected_lam_of(sig_i);
    j_hs = selected_lam_hs(sig_i);
    
    lambda1 = lambdaVals(j_s);
    lambda2 = lambdaOFVals(j_of);
    lambda3 = lambdaHSVals(j_hs);
    HSiters = 100;

    outFile = fullfile(topDir2,...
        ['results_trial_1_sig_',num2str(sig_i)],...
        ['output_j',num2str(j_s),'_',num2str(j_of),'_',num2str(j_hs),'*.mat']);
    dFiles = dir(outFile);
    loadFile = fullfile(dFiles(1).folder,dFiles(1).name);
    load(loadFile)

    [outVector, outLabels] = compute_recov_metrics(outputs,dataset,sigmas(sig_i),...
        lambda1,lambda2,lambda3,HSiters);
    outTable(:,sig_i,2) = outVector;

end

%% Plot Results
% Compute SNR
[meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);

figure
for i = 1:8
    subplot(4,2,i)
    plot(meanSNR,outTable(i,:,1),'o-')
    hold on
    plot(meanSNR,outTable(i,:,2),'s-')
    ylabel(outLabels{i})
    xlabel('SNR')
end