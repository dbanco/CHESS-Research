%% Wasserstein regularization test script

% Compute distance matrix
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
P.params.wLam = 10;
P.params.gamma = 1;
P.params.beta = 1.1;
P.params.maxIter = 1000;
P.params.maxIterReg = 200;
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

% Construct distance matrix
N = P.num_var_t*P.num_var_r;
THRESHOLD = 32;

D = ones(N,N).*THRESHOLD;
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1])
            for jj = max([1 j-THRESHOLD+1]):min([P.num_var_r j+THRESHOLD-1])
                ind1 = i + (j-1)*P.num_var_t;
                ind2 = ii + (jj-1)*P.num_var_t;
                D(ind1,ind2)= sqrt((i-ii)^2+(j-jj)^2); 
            end
        end
    end
end
D = D./max(D(:));
%% Run FISTA updating solution and error array
vdfs = load_neighbors_vdf(output_dir,baseFileName,P);
[x_hat, err_new, ~, ~] = space_wasserstein_FISTA_Circulant(A0ft_stack,b,vdfs,D,P.var_theta,P.var_rad,fileData.x_hat,P.params);
err = [fileData.err(:);err_new(:)];


%% Plot vdfs

vdf_wass = squeeze(sum(sum(fileData.x_hat,1),2));
vdf_wass = vdf_wass/sum(vdf_wass(:));

figure
subplot(3,1,1)
imagesc(vdfs{1})
title('vdf t-1')
n_vdf = vdfs{1};
dist_tm1 = sinkhornKnoppTransport(D, P.params.wLam, vdf_wass, n_vdf(:))

subplot(3,1,2)
imagesc(vdf_wass)
title('vdf W reg')

subplot(3,1,3)
imagesc(vdfs{2})
title('vdf t+1')
n_vdf = vdfs{2};
dist_tp1 = sinkhornKnoppTransport(D, P.params.wLam, vdf_wass, n_vdf(:))

vdf1 = vdfs{1};
vdf2 = vdfs{2};
dist_tm1_tp1 = sinkhornKnoppTransport(D, P.params.wLam, vdf1, vdf2(:))

[awmv_az_tm1, awmv_rad_tm1] = computeAWMV_vdf(vdf1,P.var_theta,P.var_rad)
[awmv_az_t, awmv_rad_t] = computeAWMV_vdf(vdf_wass,P.var_theta,P.var_rad)
[awmv_az_tp1, awmv_rad_tp1] = computeAWMV_vdf(vdf2,P.var_theta,P.var_rad)

%% Save outputs, updating the coefficients of the previous iteration
% save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P')
