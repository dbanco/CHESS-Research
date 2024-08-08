% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar_reduced');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','al7075_311_polar_fit_spatial_1');
mkdir(outputdir)

% Function
funcName1 = 'init_FISTA_Circulant';
funcName2 = 'wrap_space_FISTA_Circulant';
jobDir1 = fullfile(datadir,'job_al7075_311_spatial_init');
jobDir2 = fullfile(datadir,'job_al7075_311_spatial');
mkdir(jobDir1)
mkdir(jobDir2)

%% Fixed Parameters
% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;
P.sampleDims = [37,5];

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  2,    P.num_var_r).^2;

% basis weighting
P.weight = 1;
P.betap = 1;
P.alphap = 5e-2;

% fista params
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 20;
params.gamma = 1;
params.maxCycles = 10;
params.isNonnegative = 1;
P.params = params;

%% Parameters to vary
k = 0;
for img = 0:184
    for load_step = 0:4
        P.img = img;
        P.load_step = load_step;

        varin = {dataset,P,outputdir};
        funcName = funcName1;
        save(fullfile(jobDir1,['varin_',num2str(k),'.mat']),'varin','funcName')
        funcName = funcName2;
        save(fullfile(jobDir2,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

% Init script
slurm_write_bash(k-1,jobDir1,'init_batch_script.sh','0-924')
% Odd script
slurm_write_bash(k-1,jobDir2,'odd_batch_script.sh','0-924:2')
% Even script
slurm_write_bash(k-1,jobDir2,'even_batch_script.sh','1-923:2')
