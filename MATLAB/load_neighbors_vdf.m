function vdfs = load_neighbors_vdf(output_dir,baseFileName,P)
%neighbors Returns term used in gradient of vdf objective and vdfs 
%          of coefficients neighboring point (i,j)

% Get list of possible neighbors
row = ceil(P.index/(P.sampleDims(2)));
col = mod(P.index,P.sampleDims(2))+1;
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
neighbors = zeros(n,m,T,R);
vdfs = {}; 
for i = 1:size(neighbor_imgs,1)
    n_img = sub2ind(flip(P.sampleDims),neighbor_imgs(i,2),neighbor_imgs(i,1));
    n_img = n_img + P.img - P.index;
    load_file = fullfile(output_dir,sprintf(baseFileName,P.set,n_img));
    load(load_file)
    x_hat_var = x_hat;
    vdf = squeeze(sum(sum(x_hat_var)))/sum(x_hat_var(:));
    vdfs{i} = vdf;
end


