function wrap_awmv_reg_FISTA_Circulant( input_dir,P,output_dir )
%wrap_space_awmv_FISTA_Circulant Runs FISTA_Circulant loading input files and saving
% ouput files

% Load solution
baseFileName = 'fista_fit_%i_%i.mat';
fileData = load(fullfile(input_dir,sprintf(baseFileName,P.set,P.img)));
polar_image = fileData.polar_image;
P.params.t_k = fileData.P.params.t_k;
P.params.L = fileData.P.params.L;

%% Zero pad image
b = zeroPad(polar_image,P.params.zeroPad);
% Scale image by 2-norm
b = b/norm(b(:));
P.num_rad = size(b,1);
P.num_theta = size(b,2);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
    case 'norm1'
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm(P);
    case 'max'
        A0ft_stack = unshifted_basis_matrix_ft_stack(P);
end


%% Run FISTA updating solution and error array
[n_awmv_az,n_awmv_rad] = load_neighbors_awmv(input_dir,baseFileName,P);
x_init = fileData.x_hat;
[x_hat, err_new, t_k, L , ~, ~] = space_ev_FISTA_Circulant(A0ft_stack,b,n_awmv_az,n_awmv_rad,P.var_theta,P.var_rad,x_init,P.params);
err = [fileData.err(:);err_new(:)];
P.params.t_k = t_k;
P.params.L = L;

%% Save outputs, updating the coefficients of the previous iteration
save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P')


end

