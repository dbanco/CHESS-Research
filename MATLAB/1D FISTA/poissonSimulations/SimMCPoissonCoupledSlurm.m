function SimMCPoissonCoupledSlurm(snr_levels,indepPS_dir,coupledPS_dir,ps_name,mc_dir)

NN = numel(snr_levels);

% Job directory
jobDir = fullfile('/cluster','home','dbanco02',['job_',sim_name]);
mkdir(jobDir)
funcName = 'wrapCoupledPoissonMC';
k = 1;
for nn = 1:NN   
    % Make folder to save trials
    indepTrials_dir = fullfile(indepPS_dir,[mc_dir,'_',num2str(nn)]);
    coupledTrials_dir = fullfile(coupledPS_dir,[mc_dir,'_',num2str(nn)]);
    mkdir(coupledTrials_dir)

    % Load coupled selected parameter (selected_lambda2)
    ps1 = load( fullfile(coupledPS_dir,[ps_name,'_',num2str(nn)]) );
    ps2 = load(fullfile(indepPS_dir,[ps_name,'_',num2str(nn),'.mat']),'P');
    P = ps2.P;
    P.params.lambda1 = ps2.P.selected_lambdas;
    P.params.lambda2 = ps1.P.selected_lambda2;
    P.params.maxIter = 100;
    P.params.rho1 = 1.5;
    P.params.rho2 = 0.5;
    
	for i = 1:P.trials
        P.set = i;
        P.outputFile = fullfile(coupledTrials_dir,['trial_',num2str(i)]);
        varin = {P};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
	end
end
end
