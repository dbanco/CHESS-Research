%% Parameter selection
clear all
close all
disp('Setup parms')
P.set = 1;

datadir = '/cluster/shared/dbanco02/';
dataset = '/cluster/home/dbanco02/simulated_two_spot_1D_noise2_6/';
indep_dir = '/cluster/shared/dbanco02/simulated_two_spot_1D_noise2_indep_6/';
output_dir = '/cluster/shared/dbanco02/simulated_two_spot_1D_noise2_coupled_6/';

% datadir = 'D:\CHESS_data\';
% dataset = 'D:\CHESS_data\simulated_two_spot_1D_noise2_6';
% indep_dir = 'D:\CHESS_data\simulated_two_spot_1D_noise2_indep_6';
% output_dir = 'D:\CHESS_data\simulated_two_spot_1D_noise2_coupled_6';

mkdir(output_dir)
num_ims = 10;

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
load(fullfile([dataset],[prefix,'_1.mat']));
P.num_theta = size(polar_vector,1);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 20;
P.var_theta = linspace(P.dtheta/2,500,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1;
params.t_k = 1;
params.lambda = 0.08;
params.beta = 1.2;
params.maxIter = 800;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;
   
baseFileName = 'fista_fit_%i_%i.mat';

% coupled params
Pc.initialization = 'simultaneous';
Pc.preInitialized = 1;
Pc.wLam = 25;
Pc.gamma = 1;
Pc.maxIterReg = 800;
Pc.num_outer_iters = 10;
Pc.baseFileName = 'fista_fit_%i_%i.mat';
Pc.num_ims = num_ims;
Pc.prefix = 'mmpad_img';
Pc.dataset = [dataset];

% Lambda values
lambda_vals = logspace(-3,1,30); 
N = numel(lambda_vals);

% Gamma values
gamma_vals = [ 0.01,0.025 0.05,0.075,0.1,0.15,0.2]; 
M = numel(gamma_vals);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
A0 = unshifted_basis_vector_stack_norm2_zpad(P);

%% Select lambda values

disp('Selecting lambda values')

err_select = zeros(N,num_ims);
for i = 1:N
    fprintf('%i of %i \n',i,N)
    for j = 1:num_ims
        e_data = load(fullfile(indep_dir,sprintf(baseFileName,i,j)),'err');
        err_select(i,j) = e_data.err(end);
    end
end

err_select(err_select>1)=0;
imagesc(err_select)

noise_est = zeros(num_ims,1);
norm1 = zeros(num_ims,1);
norm2 = zeros(num_ims,1);
%%

for image_num = 1:num_ims
    fprintf([num2str(image_num),'\n'])
    im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
    b = squeeze(sum(im_data.polar_image,1));
    % Scale image by 2-norm
    bn = b;

    kernel = [0.0001 0.0003 0.0012 0.0042 0.0127 0.0329,...
              0.0740 0.1434 0.2399 0.3466 0.4322 0.4652];
    kernel = [kernel,fliplr(kernel(1:end-1))];
    kernel = kernel./sum(kernel(:));
    bn_hat = conv(bn,kernel,'same');

    noise_est(image_num) = norm(bn-bn_hat)/norm(bn);
    norm2(image_num) = norm(b);
    norm1(image_num) = norm(b,1);
end
% norm scaling method
noise_eta = norm2./max(norm2).*max(noise_est);

% purely residual based method
% noise_eta = noise_est;

% some other method
% noise_eta = max(noise_est)*norm2./norm1;

discrep_crit = abs(err_select'-repmat(noise_eta,1,N));
[~,lambda_indices] = min(discrep_crit');
param_select = lambda_vals(lambda_indices);
Pc.lambda_values = param_select;

%% Move independent fits to init directory

init_dir = [datadir,'mmpad_1D_coupled_simul_init'];
mkdir(init_dir)
for i = 1:num_ims
    src = fullfile(indep_dir,sprintf(baseFileName,lambda_indices(i),i));
    dest = fullfile(init_dir,sprintf(baseFileName,1,i));
    copyfile(src,dest)
end

%% Run coupled grid search

disp('Begin grid search')

for i = 1:M
    Pc.init_dir = init_dir;
    Pc.output_dirA = [datadir,'mmpad_1D_coupled_simul_',num2str(i),'a'];
    Pc.output_dirB = [datadir,'mmpad_1D_coupled_simul_',num2str(i),'b'];
    mkdir(Pc.init_dir)
    mkdir(Pc.output_dirA)
    mkdir(Pc.output_dirB)
    Pc.gamma = gamma_vals(i);
    runCoupledFISTA_1D(P,Pc)
end
