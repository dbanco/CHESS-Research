%% Parameter selection
clear all
close all
disp('Setup parms')
P.set = 1;
datadir = '/cluster/shared/dbanco02/';
dataset = ['/cluster/home/dbanco02/mmpad_polar/ring1_zero/'];
indep_dir = '/cluster/shared/dbanco02/mmpad_1D_indep_param_4/';
output_dir = '/cluster/shared/dbanco02/mmpad_1D_coupled_param_4/';

mkdir(output_dir)
num_ims = 500;


% Universal Parameters
% Ring sampling parameters
prefix = 'mmpad_img';
load([dataset,prefix,'_1.mat']);
P.num_theta = size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 15;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1;
params.t_k = 1;
params.lambda = 0.08;
params.wLam = 25;
params.beta = 1.2;
params.maxIter = 800;
params.maxIterReg = 800;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;
   
baseFileName = 'fista_fit_%i_%i.mat';

% coupled params
Pc.initialization = 'simultaneous';
Pc.preInitialized = 2;
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
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
end
A0 = unshifted_basis_vector_stack_norm2(P);

%% Select lambda values

disp('Selecting lambda values')

err_select = zeros(N,num_ims);
l0_select = zeros(N,num_ims);
l1_select = zeros(N,num_ims);
for i = 1:N
    fprintf('%i of %i \n',i,N)
    for j = 1:num_ims
        e_data = load(fullfile(indep_dir,sprintf(baseFileName,i,j)),'err','x_hat');
        err_select(i,j) = e_data.err(end-1);
        l0_select(i,j) = sum(e_data.x_hat(:) > 0);
        l1_select(i,j) = sum(e_data.x_hat(:));
    end
end
err_select(err_select > 10^10) = 0;
l0_select(l0_select > 10^10) = 0;
l1_select(l1_select > 10^10) = 0;

noise_est = zeros(num_ims,1);
norm1 = zeros(num_ims,1);
norm2 = zeros(num_ims,1);

%% Choose spartsity level and select parameters
num_coefs = 250;

sparse_crit = abs(l0_select-num_coefs);
[~,lambda_ind] = min(sparse_crit);
err_chose = zeros(num_ims,1);
l0_chose = zeros(num_ims,1);
l1_chose = zeros(num_ims,1);
lambdas_chose = zeros(num_ims,1);
inds = [lambda_ind;1:num_ims];
for i = 1:num_ims
    err_chose(i) = err_select(lambda_ind(i),i);
    l0_chose(i) = l0_select(lambda_ind(i),i);
    l1_chose(i) = l1_select(lambda_ind(i),i);
    lambdas_chose(i) = lambda_vals(lambda_ind(i));
end

figure(1)
plot(err_chose,'o')
title('error')

figure(2)
plot(l0_chose,'o')
title('l_0 norm')

figure(3)
plot(l1_chose,'o')
title('l_1 norm')

figure(4)
plot(lambdas_chose,'o')
title('\lambda values')

%% Inspect fits for chosen parameters

figure(4)
img =21;
e_data = load(fullfile(indep_dir,sprintf(baseFileName,lambda_ind(img),img)),'x_hat','polar_image');
cla
hold on
plot(squeeze(sum(e_data.polar_image,1)))
plot(Ax_ft_1D(A0ft_stack,e_data.x_hat)*norm(polar_image))
title(sprintf('err = %0.4f',err_select(lambda_ind(img),img)))
legend('data','fit')
%% Estimate noise and select parameters

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
