data_dir = '../../CHESS_data/';
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

% Plot image
figure(1)
subplot(2,1,1)
imshow(spot_org,'DisplayRange',[0 200],'Colormap',jet)
title('Original')

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


%% View basis functions used
A0_sum = squeeze(sum(sum(A0_stack_crop,1),2));
% Var coefs 
px = 6;
py = 3;

coef_sum = squeeze(sum(squeeze(sum(xhat_new,1)),1));
coef_sum_w = coef_sum.*A0_sum;
% Show spot original
figure(5)
subplot(px,py,1) 
imshow(spot_org,'DisplayRange',[0 200],'Colormap',jet)
title(['coefs sum',' --- ','weighted coefs sum' ]);
% Show spot fit
figure(5)
subplot(px,py,2) 
imshow(spot_fit,'DisplayRange',[0 200],'Colormap',jet)
title('Fit')

% Subplot coef values (fixed radial variance)
for i = 1:15
    subplot(px,py,i+2) 
    coef_slice = sum(xhat_new(:,:,i,:),4);
    imshow(coef_slice,'DisplayRange',[0 1],'Colormap',jet)
    title([num2str(sum(coef_sum(i,:))),' --- ',num2str(sum(coef_sum_w(i,:))) ]);
end
% Subplot basis function (fixed azimuthal variance)
figure(9)
for i = 1:15
    subplot(px,py,i+2) 
    basis = A0_stack_crop(:,:,i,1);
    basis_shift = shift2D(basis,10,29);
    imshow(basis_shift(1:22,1:60),'DisplayRange',[0 1],'Colormap',jet)
end

%% Run FISTA with backtracking on full image

num_theta = size(im_struc.polar_image,2);
num_rad = size(im_struc.polar_image,1);

% Generate array of 2D unshifted gaussian basis functions
A0ft_stack = unshifted_basis_matrix_ft_stack_norm(var_theta,var_rad,dtheta,drad,num_theta,num_rad);
A0_stack = unshifted_basis_matrix_stack_norm(var_theta,var_rad,dtheta,drad,num_theta,num_rad);

[xhat_new,nIter, timeSteps, errorSteps] =...
SolveFISTA_Circulant(A0ft_stack,...
                     im_struc.polar_image,...
                     'maxiteration',200,...
                     'stoppingcriterion',3,... %Objective
                     'tolerance',1e-6,...
                     'lipschitz',1e1,...
                     'lambda',50,...
                     'beta',1.2 );

params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;
[x_hat, errorSteps] = FISTA_Circulant(A0ft_stack,im_struc.polar_image,params);            
                 
                 
fit = Ax_ft_2D(A0ft_stack,xhat_new);
rel_fit_error = norm(im_struc.polar_image(:)-fit(:))/norm(im_struc.polar_image(:));
sparsity = sum(xhat_new(:)>0);

% Plot fit
figure(1)
subplot(2,1,1)
imshow(im_struc.polar_image,'DisplayRange',[0 500],'Colormap',jet)
subplot(2,1,2)
imshow(fit,'DisplayRange',[0 500],'Colormap',jet)
title('Fit')

disp(['Err ' num2str(rel_fit_error) ])
disp(['||x||_0 ' num2str(sparsity) ])

% Save test result
dir = ['fista_fit_2000_load_',...
       num2str(step), '_img_',...
       num2str(img), '.mat'];
save(dir,'xhat_new','rel_fit_error','sparsity')
% Use function to save to make this valid parfor loop
% save_vars(dir,xhat_new,rel_fit_error_new,sparse_new,...
%                        rel_fit_error_old,sparse_old);
