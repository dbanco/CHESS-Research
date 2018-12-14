function wrap_wass_reg_FISTA_Circulant( data_dir,P,output_dir )
%wrap_space_awmv_FISTA_Circulant Runs FISTA_Circulant loading input files and saving
% ouput files

% Load solution
baseFileName = 'fista_fit_%i_%i.mat';
fileData = load(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)));
polar_image = fileData.polar_image;
P = fileData.P;
%% Zero pad image
b = zeroPad(polar_image,P.params.zeroPad);
% Scale image by 2-norm
b = b/norm(b(:));
P.num_rad = size(b,1);
P.num_theta = size(b,2);

% Construct dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);

% Construct distance matrix
N = P.num_var_t*P.num_var_r;
THRESHOLD = 32;

D = ones(N,N).*THRESHOLD;
for i=1:N
    for j=max([1 i-THRESHOLD+1]):min([N i+THRESHOLD-1])
        D(i,j)= abs(i-j); 
    end
end
D = D/max(D(:)+0.1);

%% Run FISTA updating solution and error array
vdfs = load_neighbors_vdf(output_dir,baseFileName,P);
[x_hat, err_new, ~, ~] = space_wasserstein_FISTA_Circulant(A0ft_stack,b,vdfs,D,P.var_theta,P.var_rad,fileData.x_hat,P.params);
err = [fileData.err(:);err_new(:)];

%% Save outputs, updating the coefficients of the previous iteration
save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P')


end

