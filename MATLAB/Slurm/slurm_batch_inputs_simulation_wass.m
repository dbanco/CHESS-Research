% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'simulated_data_60');
prefix = 'polar_image';

% Init directory
initdir = fullfile('/cluster','shared','dbanco02','simulated_data_60_fit');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','simulated_data_60_fit_wass_reg');
mkdir(outputdir)

% Function
funcName = 'wrap_wass_reg_FISTA_Circulant';

jobDir = fullfile(datadir,['job_simulated_data_60_wass']);
mkdir(jobDir)

%% Fixed Parameters
% Ring sampling parameters
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta = size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [100,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-10;
params.L = 10;
params.lambda = 0.02;
params.wLam = 25;
params.gamma = 0.08;
params.beta = 1.2;
params.maxIter = 1000;
params.maxIterReg = 2000;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

%% Parameters to vary
img_nums = 1:100;

for img = img_nums
    P.img = img;
    P.set = 1;
    varin = {initdir,P,outputdir};
    save(fullfile(jobDir,['varin_',num2str(img),'.mat']),'varin','funcName')
end
num_jobs = 100;
slurm_write_bash(num_jobs,jobDir,'full_batch_script.sh','1-100')

