%% Load previously computed results
% results_dir = 'D:\CHESS_results\fista_fit_results\';
% load([results_dir,... 
% 'fista_out_load_',...
% num2str(step), '_img_',...
% num2str(img), '.mat']);

step = 0;
img = 35;

% Load computed results
load(['fista_fit_wavlt_',...
num2str(step),'_',...
num2str(img), '.mat']);

% Load polar_image
data_dir = 'D:\CHESS_data\al7075_311_polar\';
load([data_dir,... 
'polar_image_',...
num2str(step),'_',...
num2str(img), '.mat']);

% Compute basis matrices
A0ft_stack = unshifted_basis_matrix_ft_stack_weight(var_theta,var_rad,dtheta,drad,num_theta,num_rad,betap);
A0_stack = unshifted_basis_matrix_stack_weight(var_theta,var_rad,dtheta,drad,num_theta,num_rad,betap);

%% View full image and fit
fit_image = Ax_ft_2D(A0ft_stack,x_hat);

figure(11)
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
var_rad_idx = 5;


% Indices of region
rows = row1:row1+n-1;
cols = col1:col1+m-1;
% Subplot dimensions
px = 6;
py = 3;

% Show original and fit images restricted to region
figure(22)
subplot(px,py,1) 
imshow(log(polar_image(rows,cols)),'DisplayRange',[0 9],'Colormap',jet)
title('Original');
subplot(px,py,2) 
imshow(log(fit_image(rows,cols)),'DisplayRange',[0 9],'Colormap',jet)
title('Fit')
% Subplot coef values (fixed radial variance)
for i = 1:numel(var_theta)
    subplot(px,py,i+2) 
    coef_slice = x_hat(rows,cols,i,var_rad_idx);
    imshow(coef_slice,'DisplayRange',[0 1],'Colormap',jet)
    coef_sum = sum(coef_slice(:));
    title(num2str(coef_sum));
end

% Show original and fit images restricted to region
figure(33)
subplot(px,py,1) 
imshow(log(polar_image(rows,cols,:,:)),'DisplayRange',[0 9],'Colormap',jet)
title('Original');
subplot(px,py,2) 
imshow(log(fit_image(rows,cols,:,:)),'DisplayRange',[0 9],'Colormap',jet)
title('Fit')
% Subplot basis function 
figure(33)
for i = 1:numel(var_theta)
    subplot(px,py,i+2) 
    basis = A0_stack(:,:,i,var_rad_idx);
    basis_shift = shift2D(basis,n/2,m/2);
    imshow(basis_shift(1:n,1:m),'DisplayRange',[0 1],'Colormap',jet)
    title( ['\sigma_{\theta}/d\theta = ' num2str(sqrt(var_theta(i))/dtheta),...
            ',    \sigma_r/dr = ', num2str(sqrt(var_rad(var_rad_idx))/drad)] );
end

% Show original and fit images restricted to region
figure(44)
subplot(px,py,1) 
imshow(log(polar_image(rows,cols)),'DisplayRange',[0 9],'Colormap',jet)
title('Original');
subplot(px,py,2) 
imshow(log(fit_image(rows,cols)),'DisplayRange',[0 9],'Colormap',jet)
title('Fit')
% Subplot contributions from different coef values (fixed radial variance)
for i = 1:numel(var_theta)
    subplot(px,py,i+2) 
    signal_slice = Ax_ft_2D(A0ft_stack(:,:,i,var_rad_idx),x_hat(:,:,i,var_rad_idx));
    imshow(log(signal_slice(rows,cols)),'DisplayRange',[0 9],'Colormap',jet)
    signal_sum = sum(signal_slice(:));
    title(num2str(signal_sum));
end
