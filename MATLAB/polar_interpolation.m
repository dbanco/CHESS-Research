clear all; close all;

% Image directory
imDir = 'F:\CHESS_raw_data\al7075_mlf_mat';

% Output directory
outDir = 'F:\CHESS_data\al7075_311_polar';
mkdir(outDir)

% Ring region and sampling
center = [1023.67, 1024.61];
radius = 870;
dr = 20;
num_theta = round(2*pi*radius);
r1 = radius - dr;
r2 = radius + dr;
dtheta = 2*pi1/num_theta;

ring_params.drad = 1;
ring_params.dtheta = 1;
ring_params.num_theta = num_theta;
ring_params.num_rad = 2*dr+1;
ring_params.sampleDims = [37,5];

n_load_steps = 5;
n_imgs = 195;
for load_step = 1:n_load_steps
    for img_num = 10:n_imgs
        
        % Load full image
        load(fullfile(imDir,['al7075_mlf_',num2str(load_step),'_',num2str(img_num-1),'.mat']))
        
        % Extract ring
        polar_image = extract_ring( image,r1,r2,center,dtheta,num_theta);
        
        % Show extracted ring
        figure(1)
        imshow(log(image),'DisplayRange',[0 9],'Colormap',jet)
        hold on
        plot_circle(center,r1)
        plot_circle(center,r2)
        figure(2)
        imshow(log(polar_image),'DisplayRange',[0 9],'Colormap',jet)
        
        % Save ring
%         outFile = fullfile(outDir,['polar_image_',num2str(load_step),'_',num2str(img_num)]);
%         save(outFile,'polar_image','ring_params')
    end
end


