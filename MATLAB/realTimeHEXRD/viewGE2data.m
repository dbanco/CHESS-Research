center = [1025,1020];
r1 = 430;
r2 = 450;
% factor = 5;
% fname = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\ff_000277.ge2';
% onlineDir = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\277\';

fname = 'C:\Users\dpqb1\Documents\Data\c103_Feb\ff_001279.ge2';
onlineDir = 'C:\Users\dpqb1\Documents\Data\c103_Feb\onlineDir';
mkdir(onlineDir)
for t = 1:10
    [img_polar, rad, az] = loadGE2polar(t,fname,center,r1,r2);
    img = readGE2img(fname,t);
%     b = b./max(b(:));
%     save(fullfile(onlineDir,sprintf('polar_image_%i.mat',t)),'b')
end

imagesc(img)

for t = 1:10
    img = readGE2img(fname,t);
    save(sprintf('C:\\Users\\dpqb1\\peaknet\\data\\ge2_img%i.mat',t),'img')
end