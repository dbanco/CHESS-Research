% Read in Ge2 image
close all
T = 2;
fname = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\ff_000277.ge2';
img = readGE2img(fname,1:T);

figure
for i = 1:T
    imagesc(img(:,:,i))
    title(sprintf('image: %i',i))
end

center = [1025,1020];

% ring1
r1 = 430;
r2 = 450;

% hold on
% plot_circle(center,r1)
% plot_circle(center,r2)

num_theta = round(2*pi*r2);
dtheta = 2*pi./num_theta;
polar_image = extract_ring( img,r1,r2,center,dtheta,num_theta );

figure
imagesc(polar_image)