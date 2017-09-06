% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','al7075_311_polar_fit');


% Function
funcName = 'wrap_FISTA_Circulant';

jobDir = fullfile(datadir,'job_al7075_311');
mkdir(jobDir)

%% Parameters to vary
steps = 1:5;
imgs = 1:205;
k = 0;

for step = steps
    for img = imgs
        if ~exist(fullfile(outputdir,sprintf('fista_fit_%i_%i.mat',step,img)),'file')
            varin = {dataset,step,img,outputdir};
            save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
            k = k + 1;
        end
    end
end
        
% Job parameters
jobParams.numJobs = k;
jobParams.numWorkers = max(256,k);
jobParams.jobDir = jobDir;
jobParams.funcName = funcName;
%slurm_write_batch_script(jobParams)
% save(fullfile(jobParams.jobDir,'jobParams.mat'),'jobParams')
