% Read in Ge2 image
T = 50;
fname = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\ff_000277.ge2';
img = readGE2img(fname,1:T);

for i = 1:T
    imagesc(img(:,:,i))
    pause()
end