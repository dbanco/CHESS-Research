% Data I/O directories
data_dir = '/data/dbanco02/matlab_images/';
results_dir = '/data/dbanco02/fista_fit_results/';

% Ring sampling parameters
ring_width = 30;
num_theta= 2048;
num_rad = 2*ring_width;
dtheta = 2*pi/num_theta;
drad = 1;

% Basis function variance parameters
num_var_t = 15;
num_var_r = 10;
var_theta = linspace(dtheta,pi/32,num_var_t).^2;
var_rad   = linspace(drad,  3,       num_var_r).^2;

% Generate unshifted basis function matrices
A0ft_stack = unshifted_basis_matrix_ft_stack_norm(var_theta,var_rad,dtheta,drad,num_theta,num_rad);
A0_stack = unshifted_basis_matrix_stack_norm(var_theta,var_rad,dtheta,drad,num_theta,num_rad);

% FISTA parameters
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;


% Load polar_image
load([data_dir,... 
'polar_image_al7075_load_',...
num2str(step), '_img_',...
num2str(img), '.mat']);

% FISTA with backtracking
[x_hat, err, obj, l_0]  = FISTA_Circulant(A0ft_stack,polar_image,params);   
fit_image = Ax_ft_2D(A0ft_stack,x_hat);

% Show image and fit
figure(1)
subplot(2,1,1)
imshow(polar_image,'DisplayRange',[0 200],'Colormap',jet)
title('Original Image')
subplot(2,1,2)
imshow(fit_image,'DisplayRange',[0 200],'Colormap',jet)
title('Fit Image')


%% View basis functions used to fit

% Choose 20x60 size region to inspect
n = 30;
m = 60;
row1 = 21;
col1 = 301;
rows = row1:row1+n-1;
cols = col1:col1+m-1;

% Subplot dimensions
px = 6;
py = 3;

% Spatial sum over coefficients gets contribution of each variance basis function
coef_sum = squeeze(sum(sum(x_hat,1),2));

% Show original and fit
figure(2)
subplot(px,py,1) 
imshow(polar_image(rows,cols,:,:),'DisplayRange',[0 200],'Colormap',jet)
title(['coefs sum',' --- ','weighted coefs sum' ]);
subplot(px,py,2) 
imshow(fit_image(rows,cols,:,:),'DisplayRange',[0 200],'Colormap',jet)
title('Fit')

% Subplot coef values (fixed radial variance)
figure(2)
for i = 1:num_var_t
    subplot(px,py,i+2) 
    coef_slice = sum(xhat_new(rows,cols,i,:),4);
    imshow(coef_slice,'DisplayRange',[0 100],'Colormap',jet)
    title(num2str(sum(coef_sum(i,:))));
end
% Subplot basis function (fixed azimuthal variance)
figure(3)
for i = 1:num_var_r
    subplot(px,py,i+2) 
    basis = A0_stack(rows,cols,i,1);
    basis_shift = shift2D(basis,n/2,m/2);
    imshow(basis_shift(1:n,1:m),'DisplayRange',[0 0.1],'Colormap',jet)
    title(num2str(sum(coef_sum(:,i))));
end

