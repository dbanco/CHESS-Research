% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'simulated_data_60');
prefix = 'polar_image';

% Output directories
dirA = fullfile('/cluster','shared','dbanco02','simulated_data_60_fit_wass_reg_a');
dirB = fullfile('/cluster','shared','dbanco02','simulated_data_60_fit_wass_reg_b');
mkdir(dirA)
mkdir(dirB)

% Function
funcName1 = 'wrap_FISTA_Circulant';
funcName2 = 'wrap_wass_reg_FISTA_Circulant';

jobDir1 = fullfile(datadir,['job_simulated_data_60_wass_init']);
jobDir2 = fullfile(datadir,['job_simulated_data_60_wass_ab']);
jobDir3 = fullfile(datadir,['job_simulated_data_60_wass_ba']);
mkdir(jobDir1)
mkdir(jobDir2)
mkdir(jobDir3)

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
params.maxIter = 200;
params.maxIterReg = 5;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

%% Parameters to vary
imageDir = fullfile(dataset,prefix);

img_nums = 1:100;
P1 = P;
P1.params.tolerance = 1e-8;
k = 1;
for img = img_nums
    P1.img = img;
    P1.index = k;
    P1.set = 1;
    varin = {imageDir,P1,dirA};
    funcName = funcName1;
    save(fullfile(jobDir1,['varin_',num2str(img),'.mat']),'varin','funcName')
    k = k + 1;
end

% Regularization jobs for dirA->dirB
k = 1;
for img = img_nums
    P.img = img;
    P.index = k;
    P.set = 1;
    varin = {dirA,P,dirB};
    funcName = funcName2;
    save(fullfile(jobDir2,['varin_',num2str(k),'.mat']),'varin','funcName')
    k = k + 1;
end

% Regularization jobs for dirB->dirA
k = 1;
for img = img_nums
    P.img = img;
    P.index = k;
    P.set = 1;
    varin = {dirB,P,dirA};
    funcName = funcName2;
    save(fullfile(jobDir3,['varin_',num2str(k),'.mat']),'varin','funcName')
    k = k + 1;
end

% Init script
slurm_write_matlab(numel(img_nums),jobDir1,'parallel_simulation_wass_FISTA','batch_script.sh')
