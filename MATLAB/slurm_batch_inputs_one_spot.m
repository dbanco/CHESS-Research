% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'simulated_one_spot');
prefix = 'polar_image';

% Function
funcName = 'wrap_FISTA_Circulant';

%% Fixed Parameters
% Ring sampling parameters
load(fullfile(dataset,[prefix,'_1_1.mat']));
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [19,29];

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(0.5,20,P.num_var_t).^2;
P.var_rad   = linspace(0.5,  5,P.num_var_r).^2;
P.basis = 'norm2';

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-10;
params.L = 10;
params.lambda = 0.01;
params.beta = 1.2;
params.maxIter = 1000;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

%% Parameters to vary

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02',['one_spot_fit']);
mkdir(outputdir)

% Job directory
jobDir = fullfile('/cluster','home','dbanco02',['job_one_spot_fit']);
mkdir(jobDir)
jj = 1;
for set = 1:19
    for img = 1:29
        P.img = img;
        P.set = set;
        varin = {fullfile(dataset,prefix),P,outputdir};
        save(fullfile(jobDir,['varin_',num2str(jj),'.mat']),'varin','funcName')
        jj = jj+1;
    end
end
% slurm_write_bash(551,jobDir,'full_batch_script.sh','1-551')
slurm_write_matlab(551,jobDir,'full_batch_script.sh')
