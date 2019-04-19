function wrap_wass_weighted_reg_FISTA_Circulant( init_dir,weight_dir,P,output_dir )
%wrap_space_awmv_FISTA_Circulant Runs FISTA_Circulant loading input files and saving
% ouput files

% Load solution
baseFileName = 'fista_fit_%i_%i.mat';
fileData = load(fullfile(init_dir,sprintf(baseFileName,P.set,P.img)));
weightData = load(fullfile(weight_dir,sprintf(baseFileName,P.set,P.img)));
polar_image = fileData.polar_image;

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


% Compute data weights
wData = Ax_ft_2D(A0ft_stack,weightData.x_hat);
wData = wData/norm(wData(:));
weight_floor = prctile(wData(:),10);
weights = 1./(weight_floor + wData);

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
vdfs = load_neighbors_vdf(init_dir,baseFileName,P);
x_init = zeros(size(A0ft_stack));
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        x_init(:,:,i,j) = b/(P.num_var_t*P.num_var_r);
    end
end
[x_hat, err_new, ~, ~] = space_wasserstein_FISTA_Circulant(A0ft_stack,b,vdfs,D,x_init,P.params,weights);
err = [fileData.err(:);err_new(:)];

%% Save outputs, updating the coefficients of the previous iteration
save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P')


end

