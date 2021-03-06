%% Fixed Parameters
% dataset = 'D:\MMPAD_data';
% ringName = 'ring1_zero';
% ring_num  = 1;
% img_num = 134;
% prefix = 'mmpad_img';
% % Ring sampling parameters
% load(fullfile(dataset,ringName,[prefix,'_',num2str(img_num),'.mat']));
% P.set = ring_num;
% P.img = img_num;
% P.num_theta= size(polar_image,2);
% P.num_rad = size(polar_image,1);
% P.dtheta = 1;
% P.drad = 1;
% P.sampleDims = [546,1];
% 
% % Basis function variance parameters
% P.basis = 'norm2';
% P.num_var_t = 8;
% P.num_var_r = 12;
% P.var_theta = linspace(P.dtheta/2,6,P.num_var_t).^2;
% P.var_rad   = linspace(P.drad/2,  32,P.num_var_r).^2;
% 
% % Zero padding and mask
% maskRows = 129:133;
% zPad = [0,0];
% zMask = zeros(size(zeroPad(polar_image,zPad)));
% zMask(maskRows,:) = 1;
% zMask = onePad(zMask,zPad);
% [r,c] = find(zMask==1);
% zMask = [r,c];
% 
% % fista params
% params.stoppingCriterion = 1;
% params.tolerance = 1e-6;
% params.L = 1000;
% params.lambda = 0.1;
% params.gamma = 10;
% params.beta = 1.1;
% params.maxIter = 500;
% params.maxIterReg = 500;
% params.isNonnegative = 1;
% params.zeroPad = zPad;
% params.zeroMask = zMask;
% params.noBacktrack = 0;
% params.plotProgress = 0;
% P.params = params;

P.set = 1;
P.img = 73;

data_dir = 'D:\MMPAD_data\ring1_zero';
output_dir = 'D:\MMPAD_data\init_mmpad_reg_fit';

% Load solution
baseFileName = 'fista_fit_%i_%i.mat';
fileData = load(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)));
polar_image = fileData.polar_image;
P = fileData.P;

% fista params
P.params.stoppingCriterion = 1;
P.params.tolerance = 1e-15;
P.params.L = 500;
P.params.lambda = 0.01;
P.params.gamma = 1;
P.params.beta = 1.1;
P.params.maxIter = 1000;
P.params.maxIterReg = 1000;
P.params.isNonnegative = 1;
P.params.noBacktrack = 0;
P.params.plotProgress = 0;

%% Zero pad image
b = zeroPad(polar_image,P.params.zeroPad);
% Scale image by 2-norm
b = b/norm(b(:));
P.num_rad = size(b,1);
P.num_theta = size(b,2);

% Construct dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);

%% Run FISTA updating solution and error array
[n_awmv_az,n_awmv_rad] = load_neighbors_awmv(output_dir,baseFileName,P);
[x_hat, err_new, ~, ~] = space_ev_FISTA_Circulant(A0ft_stack,b,n_awmv_az,n_awmv_rad,P.var_theta,P.var_rad,fileData.x_hat,P.params);
err = [fileData.err(:);err_new(:)];
