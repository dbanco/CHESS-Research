%% Parameter selection
clear all
close all
P.set = 1;

num_ims = 20;
dset_name = 'gnoise4_nonorm';

% out_dir = '/cluster/shared/dbanco02/ADMM_Sherman_indep1/';
out_dir = 'D:\CHESS_data\ADMM_CG_indep1\';
mkdir(out_dir)

for jjj = 3
% dataset = ['/cluster/home/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(jjj)];
% output_dir = [out_dir,'simulated_two_spot_1D_',dset_name,'_',num2str(jjj),'_indep1/'];
dataset = ['D:\CHESS_data\simulated_two_spot_1D_',dset_name,'_',num2str(jjj)];
output_dir = [out_dir,'simulated_two_spot_1D_',dset_name,'_',num2str(jjj),'_indep\'];
mkdir(output_dir)

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta = size(polar_vector,1);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 20;
P.var_theta = linspace(P.dtheta/2,50,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% lasso/admm params
params.lambda1 = 0.01;
params.rho1 = 1;

params.tau = 1.05;
params.mu = 2;
params.adaptRho = 1;
params.alpha = 1.8;
params.maxIter = 800;
params.conjGradIter = 50;

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.tolerance = 1e-6;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.plotProgress = 0;
params.verbose = 1;

P.params = params;
   
P.baseFileName = 'fista_fit_%i.mat';

% Lambda values
lambda_vals = logspace(-4,1,30);
NumLambdas = numel(lambda_vals);
P.lambda_values = lambda_vals;

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
T = num_ims;
[N,M] = size(A0ft_stack);
%% Run grid search

% Load data
B = zeros(N,T);
for image_num = 1:num_ims
    im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
    % Zero pad image
    B(:,image_num) = im_data.polar_vector;
end

% Initialize solution
X_init = zeros([N,M,T]);

% Search lambda values
for ii = 1:NumLambdas
    P_local = P;
    P_local.params.lambda1 = lambda_vals(ii);
    P_local.set = ii;
    [X_hat, err, obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_1D(A0ft_stack,B,X_init,P_local.params);
    try
        save_output(output_dir,X_hat,err,obj,im_data.polar_image,P_local);
    catch
        save_output(output_dir,X_hat,err,obj,im_data.polar_vector,P_local);
    end
end


end

function save_output(output_dir,X_hat,err,obj,polar_image,P)
    save(fullfile(output_dir,sprintf(P.baseFileName,P.set)),'X_hat','err','obj','polar_image','P');
end

