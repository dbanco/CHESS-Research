function wrap_wass_local_reg_FISTA_Circulant( init_dir,P,output_dir )
%wrap_space_awmv_FISTA_Circulant Runs FISTA_Circulant loading input files and saving
% ouput files

% Load solution
baseFileName = 'fista_fit_%i_%i.mat';
fileData = load(fullfile(init_dir,sprintf(baseFileName,P.set,P.img)));
polar_image = fileData.polar_image;

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
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1])
            for jj = max([1 j-THRESHOLD+1]):min([P.num_var_r j+THRESHOLD-1])
                ind1 = i + (j-1)*P.num_var_t;
                ind2 = ii + (jj-1)*P.num_var_t;
                D(ind1,ind2)= sqrt((i-ii)^2+(j-jj)^2); 
            end
        end
    end
end
D = D./max(D(:));


%% Run FISTA updating solution and error array
x_init = zeros(size(A0ft_stack));
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        x_init(:,:,i,j) = b/(P.num_var_t*P.num_var_r);
    end
end
[x_hat, err_new, ~, ~] = local_wasserstein_FISTA_Circulant(A0ft_stack,b,D,x_init,P.params);
err = [fileData.err(:);err_new(:)];

%% Save outputs, updating the coefficients of the previous iteration
save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P')


end

