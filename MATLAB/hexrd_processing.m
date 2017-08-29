clear all; close all;

outDir = 'D:\CHESS_data\al7075_311_polar';
mkdir(outDir)

%% Load HEXRD image 
imDir = 'D:\CHESS_raw_data\al7075_mlf_mat';

%% AL7075 311 ring (there is also a 440 ring)
radius = 718;
dr = 20;
num_theta = 2048;

center = [1020.67, 1024.61];
r1 = radius - dr;
r2 = radius + dr;
dtheta = 2*pi/num_theta;

% figure(1)
% imshow(image,'DisplayRange',[0 200],'Colormap',jet)
% hold on
% plot_circle(center,radius)
% figure(2)
% imshow(polar_image,'DisplayRange',[0 200],'Colormap',jet)

for i = 19:205
    for j = 1:5
        img_num = i-1;
        load_step = j-1;
        % Load full image
        load(fullfile(imDir,['al7075_mlf_',num2str(load_step),'_',num2str(img_num),'.mat']))

        % Extract ring
        polar_image = extract_ring( image,r1,r2,center,dtheta,num_theta);
        
        % Save ring
        outFile = fullfile(outDir,['polar_image_',num2str(load_step),'_',num2str(img_num)]);
        save(outFile,'polar_image')
    end
end


