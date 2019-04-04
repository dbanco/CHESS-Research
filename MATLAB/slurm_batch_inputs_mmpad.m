% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'mmpad_polar');
ringName = 'ring1_zero';
ring_num  = 1;
prefix = 'mmpad_img';

% Function
funcName = 'wrap_FISTA_Circulant';

%% Fixed Parameters
% Ring sampling parameters
load(fullfile(dataset,ringName,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [546,1];

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 8;
P.var_theta = linspace(P.dtheta/2,32,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  6,P.num_var_r).^2;
P.basis = 'norm2';

% Zero padding and mask
maskCols = 129:133;
zPad = [0,0];
zMask = zeros(size(zeroPad(polar_image,zPad)));
zMask(:,maskCols) = 1;
zMask = onePad(zMask,zPad);
[r,c] = find(zMask==1);
zMask = [r,c];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-7;
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
img_nums = 1:546;

for ring_num = 1
    k = 1;
    ringName = sprintf('ring%i_zero',ring_num);
    
    % Output directory
    outputdir = fullfile('/cluster','shared','dbanco02',['mmpad_',ringName,'_fit']);
    mkdir(outputdir)

    % Job directory
    jobDir = fullfile('/cluster','home','dbanco02',['job_mmpad_',ringName,'_fit']);
    mkdir(jobDir)

    for img = img_nums
        P.img = img;
        P.set = ring_num;
        varin = {fullfile(dataset,ringName,prefix),P,outputdir};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
    slurm_write_bash(k-1,jobDir,'full_batch_script.sh','1-546')
end
