%% Parameter selection
clear all
close all
P.set = 1;

% dataset = 'D:\CHESS_data\simulated_two_spot_1D_noise2_6';
% output_dir = 'D:\CHESS_data\simulated_two_spot_1D_noise2_indep_6';

num_ims = 20;
dset_name = 'gnoise4_nonorm';
% out_dir = '/cluster/shared/dbanco02/ADMM_Sherman_indep1/';
out_dir = 'D:\CHESS_data\ADMM_Sherman_indep3\';
mkdir(out_dir)
for jjj = 3
% dataset = ['/cluster/home/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(jjj)];
% output_dir = [out_dir,'simulated_two_spot_1D_',dset_name,'_',num2str(jjj),'_indep1/'];
dataset = ['D:\CHESS_data\simulated_two_spot_1D_',dset_name,'_',num2str(jjj)];
output_dir = [out_dir,'simulated_two_spot_1D_',dset_name,'_',num2str(jjj),'_indep1\'];
mkdir(output_dir)

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
P.var_theta = linspace(P.dtheta/2,50,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% ism admm params
params.stoppingCriterion = 2;
params.lambda = 0.08;
params.rho = 4;
params.tau = 1.05;
params.mu = 10;
params.adaptRho = 1;
params.alpha = 1.8;
params.isNonnegative = 1;
params.maxIter = 1000;
params.tolerance = 1e-8;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.plotProgress = 0;
P.params = params;
   
baseFileName = 'fista_fit_%i_%i.mat';

% Lambda values
lambda_vals = logspace(-4,1,30);
N = numel(lambda_vals);
P.lambda_values = lambda_vals;

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
    bn = im_data.polar_vector;
      
    % Initial solution
    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        x_init(:,i) = bn/P.num_var_t;
    end
    
    for ii = 1:N
        P_local = P;
        P_local.params.lambda = lambda_vals(ii);
        P_local.set = ii;
        [x_hat,err,obj] = convADMM_LASSO_Sherman_1D(A0ft_stack,bn,x_init,P_local.params);
        try
            save_output(output_dir,baseFileName,x_hat,err,im_data.polar_image,P_local,image_num);
        catch
            save_output(output_dir,baseFileName,x_hat,err,im_data.polar_vector,P_local,image_num);
        end
    end
end

end

function save_output(output_dir,baseFileName,x_hat,err,polar_image,P,image_num)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','polar_image','P');
end
function save_obj(output_dir,pass,image_num,obj)
    save(fullfile(output_dir,sprintf('objective_%i_%i.mat',pass,image_num)),'obj');
end
