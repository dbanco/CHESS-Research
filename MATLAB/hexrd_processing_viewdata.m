clear all; close all;
imDir = 'D:\CHESS_raw_data\al7075_mlf_mat';

%% AL7075 311 ring
outDir = 'F:\CHESS_data\al7075_311_polar';
mkdir(outDir)
radius = 712;
dr = 20;
%% AL7075 222 ring
% outDir = 'F:\CHESS_data\al7075_222_polar';
% mkdir(outDir)
% radius = 748;
% dr = 18;
%% AL7075 400 ring
% outDir = 'F:\CHESS_data\al7075_400_polar';
% mkdir(outDir)
% radius = 870;
% dr = 20;
%% Extract polar ring images
num_theta = 2048;
center = [1023.67, 1024.61];
r1 = radius - dr;
r2 = radius + dr;
dtheta = 2*pi/num_theta;
 
for j = [1,5]
    load_step = j-1
    % Load dark image
    dark_image = load(fullfile(imDir,['al7075_mlf_',num2str(load_step),'_0.mat']));
    for i = 176
        img_num = i-1;
        % Load full image
        load(fullfile(imDir,['al7075_mlf_',num2str(load_step),'_',num2str(img_num),'.mat']))
        % Subtract dark image
        img = image - dark_image.image;
        % Clip negative values
        img(img<0) = 0;
        % Extract ring
        polar_image = extract_ring( img,r1,r2,center,dtheta,num_theta);
        
        % Plot ring (uncomment to adjust ring radius values)
        figure(1)
        imshow(log(img),'DisplayRange',[0 9],'Colormap',jet)
        hold on
        plot_circle(center,r1)
        plot_circle(center,r2)
        figure(2)
        imshow(log(polar_image),'DisplayRange',[0 9],'Colormap',jet)
        
        % Save ring
%         outFile = fullfile(outDir,['polar_image_',num2str(load_step),'_',num2str(img_num)]);
%         save(outFile,'polar_image')
    end
end


