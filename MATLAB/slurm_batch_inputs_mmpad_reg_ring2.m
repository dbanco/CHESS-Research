% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'mmpad_polar');
ringName = 'ring2_zero';
ring_num  = 2;
prefix = 'mmpad_img';

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02',['mmpad_',ringName,'_fit_reg4']);
mkdir(outputdir)

% Function
funcName1 = 'wrap_FISTA_Circulant';
funcName2 = 'wrap_awmv_reg_FISTA_Circulant';
jobDir1 = fullfile(datadir,['job_mmpad_',ringName,'_reg_init']);
jobDir2 = fullfile(datadir,['job_mmpad_',ringName,'_reg']);
mkdir(jobDir1)
mkdir(jobDir2)

%% Fixed Parameters

% Ring sampling parameters
load(fullfile(dataset,ringName,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [546,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 8;
P.num_var_r = 12;
P.var_theta = linspace(P.dtheta/2,6,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  32,P.num_var_r).^2;

% Zero padding and mask
maskRows = 129:133;
zPad = [0,0];
zMask = zeros(size(zeroPad(polar_image,zPad)));
zMask(maskRows,:) = 1;
zMask = onePad(zMask,zPad);
[r,c] = find(zMask==1);
zMask = [r,c];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-6;
params.L = 500;
params.lambda = 0.01;
params.gamma = 1;
params.beta = 1.1;
params.maxIter = 500;
params.maxIterReg = 100;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

%% Parameters to vary
img_nums = 1:546;
k = 0;
for img = img_nums
    P.img = img;
    P.set = ring_num;
    varin = {fullfile(dataset,ringName,prefix),P,outputdir};
    funcName = funcName1;
    save(fullfile(jobDir1,['varin_',num2str(k),'.mat']),'varin','funcName')
    funcName = funcName2;
    save(fullfile(jobDir2,['varin_',num2str(k),'.mat']),'varin','funcName')
    k = k + 1;
end

% Init script
slurm_write_bash(k-1,jobDir1,'init_batch_script.sh','0-545')
% Odd script
slurm_write_bash(k-1,jobDir2,'odd_batch_script.sh','0-545:2')
% Even script
slurm_write_bash(k-1,jobDir2,'even_batch_script.sh','1-545:2')
