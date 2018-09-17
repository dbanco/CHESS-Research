function wrap_awmv_reg_FISTA_Circulant( data_dir,P,output_dir )
%wrap_space_awmv_FISTA_Circulant Runs FISTA_Circulant loading input files and saving
% ouput files

% Load solution
baseFileName = 'fista_fit_%i_%i.mat';
load(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)))

%% Zero pad image
b = zeroPad(polar_image,P.params.zeroPad);
% Scale image by 2-norm
b = b/norm(b(:));
P.num_rad = size(b,1);
P.num_theta = size(b,2);

% Construct dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);

%% Run FISTA updating solution and error array
[n_awmv_az,n_awmv_rad] = load_neighbors_awmv(output_dir,baseFileName,P);
[x_hat, err_new, ~, ~] = space_ev_FISTA_Circulant(A0ft_stack,b,n_awmv_az,P.var_theta,x_hat,P.params);
err = [err(:);err_new(:)];

%% Save outputs, updating the coefficients of the previous iteration
save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P')


end

