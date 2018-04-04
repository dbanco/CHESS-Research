% Load in fit data
output_dir = 'E:\\CHESS_data\\al7075_311_polar_fit_spatial_full5\\';
% output_dir = 'E:\\CHESS_data\\al7075_311_polar_fit3\\';
baseFileName = 'spatial_fit_%i_%i.mat';
% baseFileName = 'fista_fit_%i_%i.mat';
load_step = 0;
img = 109;

% Compute vdf objective
load(fullfile(output_dir,sprintf(baseFileName,load_step,img)));
[x_neighbors,vdfs] = load_neighbors_vdf([output_dir,baseFileName],P);
x_vdf = squeeze(sum(sum(x_hat,1),2))/sum(x_hat(:));
vdf_obj = 0;
for idx = 1:numel(vdfs)
    vdf_obj = vdf_obj + norm(vdfs{idx} - x_vdf)^2;
end

vdf_obj