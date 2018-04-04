% Load in fit data
output_dir = 'E:\\CHESS_data\\al7075_311_polar_fit_spatial_sub0\\';
% output_dir = 'E:\\CHESS_data\\al7075_311_polar_fit3\\';
baseFileName = 'spatial_fit_%i_%i.mat';
% baseFileName = 'fista_fit_%i_%i.mat';
load_step = 0;
A0ft_stack = unshifted_basis_matrix_ft_stack(P);

for img = 10:20

    % Compute vdf objective
    load(fullfile(output_dir,sprintf(baseFileName,load_step,img)));
    [x_neighbors,vdfs] = load_neighbors_vdf([output_dir,baseFileName],P);
    x_vdf = squeeze(sum(sum(x_hat,1),2))/sum(x_hat(:));
    vdf_obj = 0;
    for idx = 1:numel(vdfs)
        vdf_obj = vdf_obj + norm(vdfs{idx} - x_vdf)^2;
    end
    
    % Compute normal objective
    f = 0.5*norm(polar_image-Ax_ft_2D(A0ft_stack,x_hat))^2;
    f2 = P.params.lambda * norm(x_hat(:),1); 
    
    
    est = sum(polar_image(:));

    [f f2 vdf_obj est vdf_obj*est*10000]
end