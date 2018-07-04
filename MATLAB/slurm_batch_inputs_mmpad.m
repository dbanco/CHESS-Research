% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'mmpad_polar');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','mmpad_ring1_fit');

mkdir(outputdir)
% Function
funcName = 'wrap_mmpad_norm2_FISTA_Circulant';

jobDir = fullfile('/cluster','home','dbanco02','mmpad_ring1_fit');
mkdir(jobDir)

%% Fixed Parameters

% Ring sampling parameters
P.num_theta= 51;
P.num_rad = 256;
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [546,1];

% Basis function variance parameters
P.num_var_t = 10;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,16,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  64,P.num_var_r).^2;

% fista params
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 100;
params.beta = 1.2;
params.maxIter = 300;
params.isNonnegative = 1;
P.params = params;

%% Parameters to vary
img_nums = 1:546;
ring_nums = 1:4;

k = 0;
for ring_num = ring_nums
    for img = img_nums
        P.img = img;
        P.ring_num = ring_num;
        ringName = sprintf('ring%i',ring_num);
        varin = {fullfile(dataset,ringName),P,outputdir};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

slurm_write_bash(k-1,jobDir,'full_batch_script.sh','0-2184')
