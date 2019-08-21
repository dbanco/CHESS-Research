P.set = 1;

dataset = 'D:\MMPAD_data\ring1_zero_subset';
num_ims = 10;

%% Universal Parameters
% Ring sampling parameters
prefix = 'mmpad_img';
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 15;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

%% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1000;
params.t_k = 1;
params.lambda = 0.01;
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

output_dir = '/cluster/shared/dbanco02/mmpad_1D_subset_independent';
mkdir(output_dir)

% iterate over each image
parpool(4)
parfor image_num = 1:num_ims
    im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
    %% Zero pad image
    b = zeroPad(im_data.polar_image,P.params.zeroPad);
    % Scale image by 2-norm
    b = b/norm(b(:));
    % Sum radially
    b = squeeze(sum(b,1));
    % Construct dictionary
    switch P.basis
        case 'norm2'
            A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
    end

    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
        x_init(:,i) = b/P.num_var_t;
    end

    [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D(A0ft_stack,b,x_init,P.params);

    % Output data
    save_output(output_dir,baseFileName,x_hat,err,im_data.polar_image,P,image_num);
    save_obj(output_dir,0,image_num,obj);
end

function save_output(output_dir,baseFileName,x_hat,err,polar_image,P,image_num)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','polar_image','P');
end
function save_obj(output_dir,pass,image_num,obj)
    save(fullfile(output_dir,sprintf('objective_%i_%i.mat',pass,image_num)),'obj');
end
