% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','al7075_311_polar_fit_alphap');
mkdir(outputdir)
% Function
funcName = 'wrap_FISTA_Circulant';

jobDir = fullfile(datadir,'job_al7075_311_alphap');
mkdir(jobDir)

%% Fixed Parameters
P.img = 35;
P.load_step = 0;

% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*ring_width+1;
P.dtheta = 2*pi/num_theta;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(dtheta,pi/64,num_var_t).^2;
P.var_rad   = linspace(drad,  2,       num_var_r).^2;

% basis weighting
P.weight = 1;
P.betap = P.dtheta*P.drad;

% fista params
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;
P.params = params;

%% Parameters to vary
alphaps = [1 5 10 20 50 100 200 500 1000];
k = 0;
for alphap = alphaps
    P.alphap = alphap;

    varin = {dataset,P,outputdir};
    save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
    k = k + 1;
end

slurm_write_bash(k-1,jobDir)
