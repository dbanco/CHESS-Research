% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'simulated_data_60');
prefix = 'polar_image';

% Function
funcName = 'wrap_FISTA_Circulant';

%% Fixed Parameters
% Ring sampling parameters
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [100,1];

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r).^2;
P.basis = 'norm2';

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-10;
params.L = 10;
params.lambda = 0.005;
params.beta = 1.2;
params.maxIter = 1000;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

%% Parameters to vary
img_nums = 1:100;

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02',['simulated_data2_fit']);
mkdir(outputdir)

% Job directory
jobDir = fullfile('/cluster','home','dbanco02',['job_simulated_data2_fit']);
mkdir(jobDir)

for img = img_nums
    P.img = img;
    P.set = 1;
    varin = {fullfile(dataset,prefix),P,outputdir};
    save(fullfile(jobDir,['varin_',num2str(img),'.mat']),'varin','funcName')
end
slurm_write_bash(100,jobDir,'full_batch_script.sh','1-100')

