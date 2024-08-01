% Local wrap test script
data_dir = 'D:\MMPAD_data\ring1_zero\mmpad_img';
output_dir = 'D:\MMPAD_data\test_output';
ringName = 'ring1_zero';
load([data_dir,'_30.mat']);

% Ring sampling parameters
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [546,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 8;
P.num_var_r = 12;
P.var_theta = linspace(P.dtheta/2,10,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  32,P.num_var_r).^2;

% Zero padding and mask
maskRows = 129:133;
zPad = [0,0];

zMask = zeros(size(zeroPad(polar_image,zPad)));
zMask(maskRows,:) = 1;
zMask = onePad(zMask,zPad);
[r,c] = find(zMask==1);
zMask = [r,c];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-6;
params.L = 1;
params.lambda = 0.1;
params.beta = 1.1;
params.maxIter = 500;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

P.img = 1;
P.set = 30;

%% Fit
wrap_FISTA_Circulant(data_dir,P,output_dir)

%% View Result
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
img_fit = Ax_ft_2D(A0ft_stack,x_hat);
img_fit_z = forceMaskToZero(img_fit,zMask);
rel_err = norm(polar_image(:)-img_fit_z(:))/norm(polar_image(:))
l0 = sum(x_hat(:)>0)
lim1 = 0;
lim2 = max(polar_image(:));
% Plot both images
figure(6)
subplot(1,3,1)
imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
subplot(1,3,2)
imshow(img_fit*norm(polar_image(:)),'DisplayRange',[lim1 lim2],'Colormap',jet);
subplot(1,3,3)
imshow(abs(polar_image-img_fit*norm(polar_image(:))),'DisplayRange',[lim1 2*lim2],'Colormap',jet);