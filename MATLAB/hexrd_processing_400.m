clear all; close all;
outDir = 'F:\CHESS_data\al7075_400_polar';
mkdir(outDir)

%% Load HEXRD image 
imDir = 'D:\CHESS_raw_data\al7075_mlf_mat';

%% AL7075 311 ring 
radius = 800;
dr = 20;
num_theta = 2048;

center = [1023.67, 1024.61];
r1 = radius - dr;
r2 = radius + dr;
dtheta = 2*pi/num_theta;

for i = 116:205
    for j = 1:5
        img_num = i-1;
        load_step = j-1;
        % Load full image
        load(fullfile(imDir,['al7075_mlf_',num2str(load_step),'_',num2str(img_num),'.mat']))

        % Extract ring
        polar_image = extract_ring( image,r1,r2,center,dtheta,num_theta);
        
        % Plot ring
        figure(1)
        imshow(image,'DisplayRange',[0 200],'Colormap',jet)
        hold on
        plot_circle(center,r1)
        plot_circle(center,r2)
        figure(2)
        imshow(polar_image,'DisplayRange',[0 200],'Colormap',jet)
        
        % Save ring
%         outFile = fullfile(outDir,['polar_image_',num2str(load_step),'_',num2str(img_num)]);
%         save(outFile,'polar_image')
    end
end


