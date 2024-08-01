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
P.params.tolerance = 1e-8;
P.params.L = 1;
P.params.lambda = 0.02;
P.params.gamma = 0.05;
P.params.wLam = 25;
P.params.beta = 1.5;
P.params.maxIter = 1000;
P.params.maxIterReg = 5000;
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
x_init = zeros(size(A0ft_stack));
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        x_init(:,:,i,j) = b/(P.num_var_t*P.num_var_r);
    end
end
% x_init = fileData.x_hat/10;
lambda_reg = 0.05;
P.params.lambda = lambda_reg;
[x_reg, err_new, ~, ~] = space_wasserstein_FISTA_Circulant(A0ft_stack,b,vdfs,D,x_init,P.params);


lambda_unreg = 0.05;
P.params.lambda = lambda_unreg;
[x_unreg, err_unreg, ~, ~] = FISTA_Circulant(A0ft_stack,b,x_init,P.params);

%% Display transport distances
% Definitely reduces Wass dist, but effect is weird. Need to run on whole
% dataset to avoid. Noise is a little heavy, might need to add larger basis
% functions since it wants to use largest basis functions right now. Not to
% exceed image size though.

% vdf_unreg = squeeze(sum(sum(fileData.x_hat,1),2));
vdf_unreg = squeeze(sum(sum(x_unreg,1),2));
vdf_unreg = vdf_unreg/sum(vdf_unreg(:));

vdf_wass = squeeze(sum(sum(x_reg,1),2));
vdf_wass = vdf_wass/sum(vdf_wass(:));

n_vdf1 = vdfs{1};
n_vdf2 = vdfs{2};

dist_tm1 = sinkhornKnoppTransport( vdf_wass(:), n_vdf1(:), P.params.wLam,D)
dist_tp1 = sinkhornKnoppTransport( vdf_wass(:), n_vdf2(:), P.params.wLam,D)

dist_unreg_tm1 = sinkhornKnoppTransport( vdf_unreg(:), n_vdf1(:), P.params.wLam, D)
dist_unreg_tp1 = sinkhornKnoppTransport( vdf_unreg(:), n_vdf2(:), P.params.wLam,D)
dist_tm1_tp1 = sinkhornKnoppTransport(n_vdf1(:), n_vdf2(:), P.params.wLam, D)

% Plot vdfs
figure(1)
subplot(4,1,1)
imagesc(n_vdf1)
title(['vdf unreg t-1: ',sprintf(' Wd_{t-1,t+1}= %0.4f',dist_tm1_tp1)])

subplot(4,1,2)
imagesc(vdf_wass)
title(['vdf W-reg t: ',sprintf('Wd_{t-1}= %0.4f, ',dist_tm1),sprintf( ' Wd_{t+1}= %0.4f',dist_tp1)])

subplot(4,1,3)
imagesc(vdf_unreg)
title(['vdf unreg t: ',sprintf('Wd_{t-1}= %0.4f, ',dist_unreg_tm1),sprintf( ' Wd_{t+1}= %0.4f',dist_unreg_tp1)])

subplot(4,1,4)
imagesc(n_vdf2)
title(['vdf unreg t+1: '])


% Display AWMVs
% Currently skewing AWMV upwards
[awmv_az_t, awmv_rad_t] = computeAWMV_vdf(vdf_wass,P.var_theta,P.var_rad)
[awmv_az_t_unreg, awmv_rad_t_unreg] = computeAWMV_vdf(vdf_unreg,P.var_theta,P.var_rad)
[awmv_az_tm1, awmv_rad_tm1] = computeAWMV_vdf(n_vdf1,P.var_theta,P.var_rad)
[awmv_az_tp1, awmv_rad_tp1] = computeAWMV_vdf(n_vdf2,P.var_theta,P.var_rad)

% Display fits
figure(2)
subplot(3,1,1)
imagesc(b)
title('data')

fit_unreg = Ax_ft_2D(A0ft_stack,x_unreg);
s_unreg = sum(x_unreg(:)>0);
subplot(3,1,2)
imagesc(fit_unreg)
title(['fit unreg err= ',...
    sprintf('%0.2f',err_unreg(end)),...
    sprintf(', coefs=%i',s_unreg),...
    ', \lambda = ',...
    sprintf(' %0.2f',lambda_unreg)])

fit_reg = Ax_ft_2D(A0ft_stack,x_reg);
s_reg = sum(x_reg(:)>0);
subplot(3,1,3)
imagesc(fit_reg)
title(['fit reg err= ',...
    sprintf('%0.2f',err_new(end-1)),...
    sprintf(', coefs=%i',s_reg),...
    ', \gamma = ',...
    sprintf('%0.2f',P.params.gamma),...
    ', \lambda = ',...
    sprintf(' %0.2f',lambda_reg)])

%% Save outputs, updating the coefficients of the previous iteration
% save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P')
