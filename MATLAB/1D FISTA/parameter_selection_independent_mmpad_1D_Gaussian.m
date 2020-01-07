%% Parameter selection
clear all
close all
P.set = 1;
dataset = ['/cluster/home/dbanco02/mmpad_polar/ring1_zero/'];
output_dir = '/cluster/shared/dbanco02/mmpad_1D_indep_param_8/';
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
P.num_var_t = 25;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,100,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r).^2;
% Zero padding and mask\
zPad = 50;
zMask = [1:50,(129:133)+50,(261+1+50):261+100];

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

% Lambda values
lambda_vals = logspace(-3,1,30);
N = numel(lambda_vals);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
A0 = unshifted_basis_vector_stack_norm2_zpad(P);

%% Run grid search
for image_num = 1:num_ims
    im_data = load(fullfile([dataset],[prefix,'_',num2str(image_num),'.mat']));
    % Zero pad image
    b = squeeze(sum(im_data.polar_image,1));
    % Scale image by 2-norm
    bn = b/norm(b(:));
    
    % Initial solution
    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        x_init(:,i) = bn/(P.num_var_t);
    end 
 
    parfor ii = 1:N
        P_new = P;
        P_new.params.lambda = lambda_vals(ii);
        P_new.set = ii;
        [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D(A0ft_stack,bn',x_init,P_new.params);
        save_output(output_dir,baseFileName,x_hat,err,obj,im_data.polar_image,P_new,image_num)
    end
end


function save_output(output_dir,baseFileName,x_hat,err,obj,polar_image,P,image_num)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','obj','polar_image','P')
end

