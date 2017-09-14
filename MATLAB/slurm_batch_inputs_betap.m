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

%% Fixed Parameters
P.img = 35;

% Ring sampling parameters
P.weight = 1;
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  2,       P.num_var_r).^2;

%% params
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;
P.params = params;

%% Parameters to vary
betaps = [0.01 0.1 1 10 20 50 100 1000];
steps = [0,4];
k = 0;

for load_step = steps
    for betap = betaps
        P.betap = betap;
        P.load_step = load_step;
        if ~exist(fullfile(outputdir,sprintf('fista_fit_%i_%i_betap_%2.2f.mat',P.load_step,P.img,P.betap)),'file')
            varin = {dataset,P,outputdir};
            save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
            k = k + 1;
        end
    end
end
