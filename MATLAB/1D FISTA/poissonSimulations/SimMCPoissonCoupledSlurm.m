function SimMCPoissonCoupledSlurm(snr_levels,indepPS_dir,coupledPS_dir,ps_name,mc_dir,sim_name)

NN = numel(snr_levels);

% Job directory
jobDir = fullfile('/cluster','home','dbanco02',['job_',sim_name,'MC']);
mkdir(jobDir)
funcName = 'wrapCoupledPoissonMC';
k = 1;
for nn = 1:NN   
    % Make folder to save trials
    indepTrials_dir = fullfile(indepPS_dir,[mc_dir,'_',num2str(nn)]);
    load(fullfile(indepTrials_dir,['trial_',num2str(1)]),'P');
    P.indepTrials_dir = indepTrials_dir;
    coupledTrials_dir = fullfile(coupledPS_dir,[mc_dir,'_',num2str(nn)]);
    mkdir(coupledTrials_dir)

    % Load coupled selected parameter (selected_lambda2)
    ps2 = load(fullfile(indepPS_dir,[ps_name,'_',num2str(nn),'.mat']),'P');
    MM = numel(P.lambda2_values);
    
    mse_all = zeros(MM,1);
    l1_norm_all = zeros(MM,1);
    tv_penalty_all = zeros(MM,1);
    awmv_err_coupled = zeros(MM,1);
    for i = 1:MM
        ps1 = load(fullfile(coupledPS_dir,...
            [ps_name,'_',num2str(nn),'_',num2str(i),'.mat']),...
            'mse','l1_norm','tv_penalty','P','X_hat','awmv_rmse');
        mse_all(i) = ps1.mse;
        l1_norm_all(i) = ps1.l1_norm;
        tv_penalty_all(i) = ps1.tv_penalty;
        awmv = computeAWMV_1D(ps1.X_hat,ps1.P.var_theta);
        awmv_err_coupled(i) = norm(ps1.P.theta_stds(:)-awmv(:))/norm(ps1.P.theta_stds(:));
    end
    [~,sel_ind] = min(awmv_err_coupled);

    
    P.params.lambda1 = ps2.P.selected_lambdas;
    P.params.lambda2 = ps2.P.lambda2_values(sel_ind);P.params.lambda2
    P.params.maxIter = 100;
    P.params.rho1 = 1;
    P.params.rho2 = 1;
    
	for i = 1:P.trials
        P.set = i;
        P.outputFile = fullfile(coupledTrials_dir,['trial_',num2str(i)]);
        varin = {P};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
	end
end
slurm_write_bash(k-1,jobDir,'full_batch_script.sh',['1-',num2str(NN*P.trials)]) %,['1-',num2str(M)])
end
