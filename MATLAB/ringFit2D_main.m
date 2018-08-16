%% Problem parameters
% Data I/O directories
data_dir = 'D:\CHESS_data\al7075_311_polar\';
results_dir = 'D:\CHESS_results\test\';

% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;

P.load_step = 4;
P.img = 35;

% Basis function dictionary parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta, pi/64, P.num_var_t).^2;
P.var_rad   = linspace(P.drad, 2, P.num_var_r).^2;

% FISTA parameters
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

% Load polar_image
load([data_dir,... 
'polar_image_',...
num2str(P.load_step),'_',...
num2str(P.img), '.mat']);

% Define zero mask
maskRows = 129:133;
zMask = zeros(size(polar_image));
zMask(maskRows,:) = 1;

% Zero pad image and update mask
zPad = [5,0];
zMask = onePad(zMask,zPad);
b = zeroPad(polar_image,P.params.zeroPad);
[r,c] = find(zMask==1);
zMask = [r,c];

% Redefine image dimensions
P.num_rad = size(b,1);
P.num_theta = size(b,2);

% Construct dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
A0_stack = unshifted_basis_matrix_stack_norm2(P);

% Initialize solution
x_init = rand(size(A0ft_stack));
x_init = x_init/norm(x_init(:));
x_init = forceMaskToZeroArray(x_init,P.zeroMask);

% Normalize image
b = b/norm(b(:));

%% FISTA with backtracking
[x_hat, err, obj, l_0]  = FISTA_Circulant(A0ft_stack,b,x_init,params);
err(end)

%% View fit
img_fit = Ax_ft_2D(A0ft_stack,x_hat);
img_fit_z = forceMaskToZero(img_fit,zMask);
lim1 = 0;
lim2 = max(polar_image(:));
% Plot data, fit, and masked fit
figure(1)
subplot(1,3,1)
imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
subplot(1,3,2)
imshow(img_fit,'DisplayRange',[lim1 lim2],'Colormap',jet);
subplot(1,3,3)
imshow(abs(polar_image-img_fit),'DisplayRange',[lim1 2*lim2],'Colormap',jet);

%% Plot error, objective, sparsity
figure(2)
subplot(3,1,1)
semilogy(err)
subplot(3,1,2)
semilogy(obj)
subplot(3,1,3)
semilogy(l_0)