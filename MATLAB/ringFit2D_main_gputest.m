data_dir = '../../CHESS_data/matlab_polar_images/';
dr = 30;

num_var_t = 15;
num_var_r = 10;

dtheta = 2*pi/2048;
drad = 1;

% Define radial and azimuthal variances
var_theta = linspace(dtheta,pi/32,num_var_t).^2;
var_rad   = linspace(1*drad,  10, num_var_r).^2;

step = 0;
img = 30;

% Load polar_image
im_struc = load([data_dir,... 
'polar_image_al7075_load_',...
num2str(step), '_img_',...
num2str(img), '.mat']);

%% Setup FISTA on cropped image
% Crop polar image
spot_org = im_struc.polar_image(20:41,1203:1262);
num_theta= size(spot_org,2);
num_rad = size(spot_org,1);


% Generate array of 2D unshifted gaussian basis functions
A0ft_stack_crop = unshifted_basis_matrix_ft_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad);
A0_stack_crop = unshifted_basis_matrix_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad);

%% Run FISTA with backtracking on cropped image
[xhat_new,nIter, timeSteps, errorSteps] =...
SolveFISTA_Circulant(A0ft_stack_crop,...
                     spot_org,...
                     'maxiteration',1000,...
                     'stoppingcriterion',3,...
                     'tolerance',1e-6,...
                     'lambda',50,...
                     'beta',1.2);

spot_fit = Ax_ft_2D(A0ft_stack_crop,xhat_new)
rel_fit_error = norm(spot_org(:)-spot_fit(:))/norm(spot_org(:))
sparsity = sum(xhat_new(:)>0)

% Plot fit
figure(1)
subplot(2,1,2)
imshow(spot_fit,'DisplayRange',[0 200],'Colormap',jet)
title('Fit')

disp(['Err ' num2str(rel_fit_error) ])
disp(['||x||_0 ' num2str(sparsity) ])

