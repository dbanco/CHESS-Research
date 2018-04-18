function wrap_space_ev_FISTA_Circulant( data_dir,P,output_dir )
%init_FISTA_Circuleant Runs FISTA_Circulant loading input files and saving
% ouput files

% Initialize solution
baseFileName = 'spatial_fit_%i_%i.mat';
load(fullfile(output_dir,sprintf(baseFileName,P.load_step,P.img)))

% Construct dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);

%% Run FISTA updating solution and error array
neighbors_ev = load_neighbors_ev_az(fullfile(output_dir,baseFileName),P);
[x_hat, err_new, ~, ~] = space_ev_FISTA_Circulant(A0ft_stack,polar_image,neighbors_ev,P.var_theta,x_hat,P.params);
err = [err(:);err_new(:)];

%% Save outputs, updating the coefficients of the previous iteration
save(fullfile(output_dir,sprintf(baseFileName,P.load_step,P.img)),'x_hat','err','polar_image','P')


end

