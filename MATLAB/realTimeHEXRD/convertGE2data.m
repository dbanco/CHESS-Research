center = [1025,1020];
r1 = 430;
r2 = 450;
% factor = 5;
% fname = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\ff_000277.ge2';
% onlineDir = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\277\';

fname = 'C:\Users\dpqb1\Documents\Data\c103_Feb\ff_000807.ge2';
onlineDir = 'C:\Users\dpqb1\Documents\Data\c103_Feb\onlineDir';
mkdir(onlineDir)
for t = 1
    [img_polar, rad, az] = loadGE2polar(t,fname,center,r1,r2);
    img = readGE2img(fname,t);
%     b = b./max(b(:));
%     save(fullfile(onlineDir,sprintf('polar_image_%i.mat',t)),'b')
end

img_rad = zeros(size(b));
img_az = zeros(size(b));
for i = 1:numel(rad)
    for j = 1:numel(az)
        img_rad(i,j) = rad(i);
        img_az(i,j) = az(j);
    end
end
true_center = [1025,1020];
detector_dist = 4261.46;


mask = zeros(2048,2048);
n = 1;
for i = 1:2048
    for j = 1:2048
        r = sqrt( (i-center(1))^2 + (j-center(2))^2 );
        theta = atan2(i-center(1),j-center(2));
        if (r1 <= r) && (r <= r2)
            mask(i,j) = 1;
            pixels(n,:) = [i,j];
            coords(n,:) = [r,theta];
            n = n + 1;
        end
    end
end