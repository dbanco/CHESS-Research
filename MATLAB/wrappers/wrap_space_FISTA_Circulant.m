function wrap_space_FISTA_Circulant( data_dir,P,output_dir )
%init_FISTA_Circuleant Runs FISTA_Circulant loading input files and saving
% ouput files

% Initialize solution
baseFileName = 'spatial_fit_%i_%i.mat';
fileData = load(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)));
x_hat = fileData.x_hat;
polar_image = fileData.polar_image;

% Construct dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack(P);

%% Run FISTA updating solution and error array
[vdfs] = load_neighbors_vdf(fullfile(output_dir,baseFileName),P);
[x_hat, err_new, ~, ~] = space_FISTA_Circulant(A0ft_stack,polar_image,x_neighbors,vdfs,x_hat,P.params);
err = [err(:);err_new(:)];

%% Save outputs, updating the coefficients of the previous iteration
save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P')


end

