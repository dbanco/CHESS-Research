P.set = 1;
P.img = 1;

dataset = '/cluster/home/dbanco02/simulated_data_two_spot_growth_far/';
output_dir = '/cluster/shared/dbanco02/two_spot_growth_far_independent';
num_ims = 25;
mkdir(output_dir)
prefix = 'polar_image';

%% Universal Parameters
% Ring sampling parameters
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r).^2;

% Zero padding and mask
zPad = [0,0];
zMask = [];

%% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1000;
params.t_k = 1;
params.lambda = 0.03;
params.wLam = 25;
params.gamma = 0.25;
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

parfor image_num = 1:num_ims
    im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
    %% Zero pad image
    b = zeroPad(im_data.polar_image,P.params.zeroPad);
    % Scale image by 2-norm
    b = b/norm(b(:));

    % Construct dictionary
    switch P.basis
        case 'norm2'
            A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
        case 'norm1'
            A0ft_stack = unshifted_basis_matrix_ft_stack_norm(P);
        case 'max'
            A0ft_stack = unshifted_basis_matrix_ft_stack(P);
    end

    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        for j = 1:P.num_var_r
            x_init(:,:,i,j) = b/(P.num_var_t*P.num_var_r);
        end
    end
  
    % Run FISTA updating solution and error array

    [x_hat,err,obj,~,~,~] = FISTA_Circulant(A0ft_stack,b,x_init,P.params);
    
    save_output(output_dir,baseFileName,x_hat,err,im_data.polar_image,P,image_num);
    save_obj(output_dir,0,image_num,obj);
end

function save_output(output_dir,baseFileName,x_hat,err,polar_image,P,image_num)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','polar_image','P');
end
function save_obj(output_dir,pass,image_num,obj)
    save(fullfile(output_dir,sprintf('objective_%i_%i.mat',pass,image_num)),'obj');
end
