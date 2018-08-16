function [n_awmv_az,n_awmv_rad] = load_neighbors_awmv(output_dir,baseFileName,P)
%neighbors Returns term used in gradient of vdf objective and vdfs 
%          of coefficients neighboring point (i,j)

% Get list of possible neighbors
row = floor(P.img/5)+1;
col = mod(P.img,5)+1;
rows = [row-1;
        row+1; 
        row;
        row];
cols = [col;
        col; 
        col-1;
        col+1];

% Remove invalid neighbors
keep = (rows >= 1).*(rows <=P.sampleDims(1)).*...
       (cols >= 1).*(cols <=P.sampleDims(2));
neighbor_imgs = [rows(keep==1),cols(keep==1)];

%Get contributions and vdfs from neighbors
n = P.num_rad;
m = P.num_theta;
T = P.num_var_t;
R = P.num_var_r;
n_awmv_az = 0;
n_awmv_rad = 0; 

for i = 1:size(neighbor_imgs,1)
    n_img = sub2ind(flip(P.sampleDims),neighbor_imgs(i,2),neighbor_imgs(i,1));
    load_file = fullfile(output_dir,sprintf(baseFileName,P.set,n_img-1));
    load(load_file)
    x_hat_var = x_hat;
    [awmv_az, awmv_rad] = computeAWMV(x_hat_var,P.var_theta,P.var_rad);
    n_awmv_az = n_awmv_az + awmv_az;
    n_awmv_rad = n_n_awmv_rad+ awmv_rad;
end
n_awmv_az = n_awmv_az/size(neighbor_imgs,1);
n_awmv_rad = n_awmv_rad/size(neighbor_imgs,1);
