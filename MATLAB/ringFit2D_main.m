%% Problem parameters
% Data I/O directories
data_dir = 'D:\CHESS_data\al7075_311_polar\';
results_dir = 'D:\CHESS_results\fista_fit_results\';

% Ring sampling parameters
ring_width = 20;
num_theta= 2048;
num_rad = 2*ring_width+1;
dtheta = 2*pi/num_theta;
drad = 1;

% Basis function variance parameters
num_var_t = 15;
num_var_r = 10;
var_theta = linspace(dtheta,pi/64,num_var_t).^2;
var_rad   = linspace(drad,  2,       num_var_r).^2;

% Generate unshifted basis function matrices
A0ft_stack = unshifted_basis_matrix_ft_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad);
A0_stack = unshifted_basis_matrix_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad);

% FISTA parameters
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;

step = 0;
img = 35;
% Load polar_image
load([data_dir,... 
'polar_image_',...
num2str(step),'_',...
num2str(img), '.mat']);

%% FISTA with backtracking
[x_hat, err, obj, l_0]  = FISTA_Circulant(A0ft_stack,polar_image,params);   

%% Or load previously computed results
load([results_dir,... 
'fista_out_load_',...
num2str(step), '_img_',...
num2str(img), '.mat']);

%% View full image and fit
fit_image = Ax_ft_2D(A0ft_stack,x_hat);

figure(1)
subplot(2,1,1)
imshow(log(polar_image),'DisplayRange',[0 9],'Colormap',jet)
title('Original Image')
subplot(2,1,2)
imshow(log(fit_image),'DisplayRange',[0 9],'Colormap',jet)
title('Fit Image')

%% Analyze regression results:  view coefficients and corresponding basis 
%  functions used to fit small region of the polar image

% Choose n x m size region to inspect
n = 20;
m = 100;
% Set starting row and col of region
row1 = 10;
col1 = 1301;
% Select radial variance of basis function (1-10)
var_rad_idx = 1;


% Indices of region
rows = row1:row1+n-1;
cols = col1:col1+m-1;
% Subplot dimensions
px = 6;
py = 3;

% Show original and fit images restricted to region
figure(2)
subplot(px,py,1) 
imshow(polar_image(rows,cols),'DisplayRange',[0 200],'Colormap',jet)
title('Original');
subplot(px,py,2) 
imshow(fit_image(rows,cols),'DisplayRange',[0 200],'Colormap',jet)
title('Fit')
% Subplot coef values (fixed radial variance)
for i = 1:num_var_t
    subplot(px,py,i+2) 
    coef_slice = x_hat(rows,cols,i,var_rad_idx);
    imshow(coef_slice,'DisplayRange',[0 1],'Colormap',jet)
    coef_sum = sum(coef_slice(:));
    title(num2str(coef_sum));
end

% Show original and fit images restricted to region
figure(3)
subplot(px,py,1) 
imshow(polar_image(rows,cols,:,:),'DisplayRange',[0 200],'Colormap',jet)
title('Original');
subplot(px,py,2) 
imshow(fit_image(rows,cols,:,:),'DisplayRange',[0 200],'Colormap',jet)
title('Fit')
% Subplot basis function 
figure(3)
for i = 1:num_var_t
    subplot(px,py,i+2) 
    basis = A0_stack(:,:,i,var_rad_idx);
    basis_shift = shift2D(basis,n/2,m/2);
    imshow(basis_shift(1:n,1:m),'DisplayRange',[0 1],'Colormap',jet)
    title( ['\sigma_{\theta}/d\theta = ' num2str(sqrt(var_theta(i))/dtheta),...
            ',    \sigma_r/dr = ', num2str(sqrt(var_rad(var_rad_idx))/drad)] );
end

