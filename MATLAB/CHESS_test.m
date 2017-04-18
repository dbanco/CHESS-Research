data_dir = 'E:\CHESS_data\matlab_polar_images\';
results_dir = 'E:\CHESS_results\matlab\';
error_dir = 'E:\CHESS_results\matlab_error\';
prune_dir = 'E:\CHESS_results\matlab_prune\';

dr = 30;

num_theta= 2048;
num_rad = 2*dr;

num_var_t = 15;
num_var_r = 10;

dtheta = 2*pi/num_theta;
drad = 1;

var_theta = linspace(dtheta,pi/32,num_var_t).^2;
var_rad   = linspace(drad,  3,       num_var_r).^2;

A0ft_stack = unshifted_basis_matrix_ft_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad);
A0_stack = unshifted_basis_matrix_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad);


%% Test original FISTA algorithm
x_hat_new3 = fista_circulant_2D(A0ft_stack, polar_image, x_hat_new, 6e5, 1, 100, eps, 1, 3, 1);
sum(x_hat_new(:)>0)/numel(x_hat_new)
fit = Ax_ft_2D(A0ft_stack,x_hat_new3);
figure(1)
imshow([fit; 400*ones(1,2048) ;polar_image],'DisplayRange',[0 400],'Colormap',jet)
%% Test FISTA with backtracking 
step = 0;
img = 20;

% Load error
err_struc = load([error_dir,... 
    'ringModel_error_load_',...
    num2str(step), '_img_',...
    num2str(img), '.mat']);

err_struc.rel_fit_error

% Load polar image
load([data_dir,... 
'polar_image_al7075_load_',...
num2str(step), '_img_',...
num2str(img), '.mat']);


% Load xhat
coef_struc = load([results_dir,... 
'ringModel_coefs_load_',...
num2str(step), '_img_',...
num2str(img), '.mat']);

[x_hat_pos2,nIter, timeSteps, errorSteps] = SolveFISTA_Circulant(A0ft_stack,...
                                            polar_image,...
                                            'maxiteration',500,...
                                            'stoppingcriterion',3,...
                                            'tolerance',1e-6,...
                                            'initialization',coef_struc.xhat);

fit = Ax_ft_2D(A0ft_stack,x_hat_pos2);
figure(1)
imshow([fit; 400*ones(1,2048) ;polar_image],'DisplayRange',[0 400],'Colormap',jet)

%% Test FISTA pruning for high error images
for step = 1:5
    parfor img = 11:200
        disp(['img ' num2str(img)])
        % Load rel_fit_error
        err_struc = load([error_dir,... 
                'ringModel_error_load_',...
                num2str(step), '_img_',...
                num2str(img), '.mat']);
        if(err_struc.rel_fit_error > 0.35)
            % Load polar_image
            im_struc = load([data_dir,... 
            'polar_image_al7075_load_',...
            num2str(step), '_img_',...
            num2str(img), '.mat']);

            % Load xhat
            coef_struc = load([results_dir,... 
            'ringModel_coefs_load_',...
            num2str(step), '_img_',...
            num2str(img), '.mat']);

            % Run 500 more iterations
            xhat_new = fista_circulant_2D(A0ft_stack,...
                                            im_struc.polar_image,...
                                            coef_struc.xhat,...
                                            6.5e5, 1, 500, eps, 1, 0, 1);
            fit_new = Ax_ft_2D(A0ft_stack,xhat_new);  
            fit_old = Ax_ft_2D(A0ft_stack,coef_struc.xhat);

            rel_fit_error_new = norm(im_struc.polar_image(:)-fit_new(:))/norm(im_struc.polar_image(:));
            rel_fit_error_old = norm(im_struc.polar_image(:)-fit_old(:))/norm(im_struc.polar_image(:));

            sparse_new = sum(xhat_new(:)>0);
            sparse_old = sum(coef_struc.xhat(:)>0);

            disp(['Err ' num2str(rel_fit_error_new) ' to ' num2str(rel_fit_error_old)])
            disp(['||x||_0 ' num2str(sparse_new) ' to ' num2str(sparse_old)])
            % Save result
            dir = [prune_dir,...
                'prune_results_load_',...
                num2str(step), '_img_',...
                num2str(img), '.mat'];
            
            save_vars(dir,xhat_new,rel_fit_error_new,sparse_new,...
                               rel_fit_error_old,sparse_old);
        end
    end
end

%% Try to improve FISTA fit
step = 0;
img = 20;

% Load error
err_struc = load([error_dir,... 
    'ringModel_error_load_',...
    num2str(step), '_img_',...
    num2str(img), '.mat']);

err_struc.rel_fit_error

% Load polar image
load([data_dir,... 
'polar_image_al7075_load_',...
num2str(step), '_img_',...
num2str(img), '.mat']);


% Load xhat
coef_struc = load([results_dir,... 
'ringModel_coefs_load_',...
num2str(step), '_img_',...
num2str(img), '.mat']);

x_hat_new1 = fista_circulant_2D(A0ft_stack, polar_image, 0, 6e5, 1, 100, eps, 1, 3, 1);
sum(x_hat_new1(:)>0)/numel(x_hat_new1)
fit = Ax_ft_2D(A0ft_stack,x_hat_new1);
figure(1)
imshow([fit; 400*ones(1,2048) ;polar_image],'DisplayRange',[0 400],'Colormap',jet)


%% Test Least Squares Initialized FISTA
step = 0;
img = 20;

% Load error
err_struc = load([error_dir,... 
    'ringModel_error_load_',...
    num2str(step), '_img_',...
    num2str(img), '.mat']);

err_struc.rel_fit_error

% Load polar image
load([data_dir,... 
'polar_image_al7075_load_',...
num2str(step), '_img_',...
num2str(img), '.mat']);

tol = 1e-8; 
maxit = 200;

A = @(x,type) Ax_Atx_ft_2D(x,type,A0ft_stack);
b = reshape(polar_image,numel(polar_image),1);
x_ls = lsqr( A,b,tol,maxit );
x_ls = reshape(x_ls,size(A0ft_stack));

fit = Ax_ft_2D(A0ft_stack,x_ls);
figure(1)
imshow([fit; 400*ones(1,2048) ;polar_image],'DisplayRange',[0 400],'Colormap',jet)
error = norm(polar_image-fit)/norm(polar_image)