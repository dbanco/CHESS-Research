function SimParamSearchCoupledPoissonSlurm(P,MM,snr_levels,...
                  output_dir,indep_dir,file_name,sim_name)
mkdir(output_dir)

% Fits for different parameters/noise levels
NN = numel(snr_levels);
P.params.maxIter = 100;
P.params.rho1 = 1;
P.params.rho2 = 1;

% Job directory
jobDir = fullfile('/cluster','home','dbanco02',['job_',sim_name]);
mkdir(jobDir)
funcName = 'wrapCoupledPoisson';
k = 1;
for nn = 1:NN
    i_name = [file_name,'_',num2str(nn),'.mat'];
    indepFile = fullfile(indep_dir,i_name);
    P.indepFile = indepFile;
    for i = 1:MM
        P.set = i;
        f_name =  [file_name,'_',num2str(nn),'_',num2str(i),'.mat'];
        P.outputFile = fullfile(output_dir,f_name);
        varin = {P};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end
slurm_write_bash(k-1,jobDir,'full_batch_script.sh',['1-',num2str(NN*MM)]) %,['1-',num2str(M)])


end
