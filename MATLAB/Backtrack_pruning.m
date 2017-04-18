data_dir = '/data/dbanco02/matlab_images/';
results_dir = '/data/dbanco02/';
error_dir = '/data/dbanco02/matlab_error/';
backtrack_dir = '/data/dbanco02/backtrack/';

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

%% Backtrack Pruning
step = 0;
parfor img = 11:200
    disp(['img ' num2str(img)])
    % Load rel_fit_error
    err_struc = load([error_dir,... 
            'ringModel_error_load_',...
            num2str(step), '_img_',...
            num2str(img), '.mat']);
    if(err_struc.rel_fit_error > 0.30)
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

        % Run 500 iterations of FISTA with backtracking
        [xhat_new,nIter, timeSteps, errorSteps] =...
        SolveFISTA_Circulant(A0ft_stack,...
                             im_struc.polar_image,...
                             'maxiteration',500,...
                             'stoppingcriterion',3,...
                             'tolerance',1e-6,...
                             'initialization',coef_struc.xhat);

        fit_new = Ax_ft_2D(A0ft_stack,xhat_new);  
        fit_old = Ax_ft_2D(A0ft_stack,coef_struc.xhat);

        rel_fit_error_new = norm(im_struc.polar_image(:)-fit_new(:))/norm(im_struc.polar_image(:));
        rel_fit_error_old = norm(im_struc.polar_image(:)-fit_old(:))/norm(im_struc.polar_image(:));

        sparse_new = sum(xhat_new(:)>0);
        sparse_old = sum(coef_struc.xhat(:)>0);

        disp(['Err ' num2str(rel_fit_error_new) ' to ' num2str(rel_fit_error_old)])
        disp(['||x||_0 ' num2str(sparse_new) ' to ' num2str(sparse_old)])
        % Save result
        dir = [backtrack_dir,...
            'prune_results_load_',...
            num2str(step), '_img_',...
            num2str(img), '.mat'];

        % Use function to save to make this valid parfor loop
        save_vars(dir,xhat_new,rel_fit_error_new,sparse_new,...
                               rel_fit_error_old,sparse_old);
    end
end
