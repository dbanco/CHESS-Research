% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'mmpad_polar');

% Function
funcName = 'wrap_mmpad_norm2_FISTA_Circulant';

%% Fixed Parameters

% Ring sampling parameters
P.num_theta= 51;
P.num_rad = 256;
P.dtheta = 0.25;
P.drad = 0.25;
P.sampleDims = [546,1];

% Basis function variance parameters
P.num_var_t = 10;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,6,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  10,P.num_var_r).^2;

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-6;
params.L = 1e3;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;
P.params = params;

%% Parameters to vary
img_nums = 1:546;


for ring_num = 1:4
    k = 0;
    ringName = sprintf('ring%i',ring_num);
    
    % Output directory
    outputdir = fullfile('/cluster','shared','dbanco02',['mmpad_',ringName,'_fit']);
    mkdir(outputdir)

    % Job directory
    jobDir = fullfile('/cluster','home','dbanco02',['job_mmpad_',ringName,'_fit']);
    mkdir(jobDir)

    for img = img_nums
        P.img = img;
        P.ring_num = ring_num;
        varin = {fullfile(dataset,ringName),P,outputdir};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
    slurm_write_bash(k-1,jobDir,'full_batch_script.sh','0-545')
end
