%% Problem parameters
% Data I/O directories
% data_dir = 'D:\CHESS_data\al7075_311_polar\';
% results_dir = 'D:\CHESS_results\test\';

data_dir = '';
results_dir = '';

% Select image
P.set = 4;
P.img = 35;

% Load polar_image
load([data_dir,... 
'polar_image_al7075_load_',...
num2str(P.set),'_img_',...
num2str(P.img), '.mat']);

% Define zero pad
zPad = [];

% Define zero mask
% maskRows = 129:133;
% zMask = zeros(size(polar_image));
% zMask(maskRows,:) = 1;
% zMask = onePad(zMask,zPad);
% [r,c] = find(zMask==1);
% zMask = [r,c];
zMask = [];

% Pad image
b = zeroPad(polar_image,zPad);


% Redefine image dimensions
P.num_rad = size(b,1);
P.num_theta = size(b,2);

%% Ring sampling parameters
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;

% Basis function dictionary parameters
P.num_var_t = 12;
P.num_var_r = 8;
P.var_theta = linspace(1, 20, P.num_var_t).^2;
P.var_rad   = linspace(1, 4, P.num_var_r).^2;

% FISTA parameters
params.stoppingCriterion = 1;
params.tolerance = 1e-6;
params.L = 10000;
params.lambda = 0.01;
params.beta = 1.1;
params.maxIter = 500;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 1;
P.params = params;

% Construct dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
A0_stack = unshifted_basis_matrix_stack_norm2(P);

% View basis functions
figure(1)
basis_sum = zeros(size(b));
shift_r = 0;
shift_c = 0;
for i = 1:P.num_var_t
    shift_r = shift_r + 5*sqrt(P.var_theta);
    for j = 1:P.num_var_r
        shift_c = shift_c + 5*sqrt(P.var_rad);
        basis_func = shift2D(squeeze(A0_stack(:,:,i,j)),round(P.num_rad/2),round(P.num_theta/2));
        up = round(P.num_theta/2+2*sqrt(P.var_theta(i)));
        dwn = round(P.num_theta/2-2*sqrt(P.var_theta(j)));
        basis_func = basis_func(:,dwn:up);
        %basis_sum = basis_sum + basis_func; 
        imshow(basis_func,'DisplayRange',[0 0.8],...
               'Colormap',jet,...
               'InitialMagnification','fit');
        pause(0.2)
    end
end


%% FISTA with backtracking

% Initialize solution
x_init = zeros(size(A0ft_stack));
x_init = forceMaskToZeroArray(x_init,params.zeroMask);

% Normalize image
% b = b/norm(b(:));
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
imshow(abs(polar_image-img_fit),'DisplayRange',[lim1 2*lim2],...
    'Colormap',jet,...
    'InitialMagnification','fit'););

%% Plot error, objective, sparsity
figure(2)
subplot(3,1,1)
semilogy(err)
subplot(3,1,2)
semilogy(obj)
subplot(3,1,3)
semilogy(l_0)