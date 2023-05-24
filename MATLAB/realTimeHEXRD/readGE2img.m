function img = readGE2img(fname,imagenum)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ge2File = fopen(fname);

% head = 'GE detector raw data';
fseek(ge2File,18,'bof');

sizexy = [2048,2048];
Npix = sizexy(1)*sizexy(2);
img = zeros([sizexy,numel(imagenum)]);

j = 1;
for i = imagenum
    pos = 4096 + (i-1)*2*Npix;
    success = fseek(ge2File,pos,'bof');
    if success == 0
        image = (fread(ge2File,Npix,"int16"));
    else
        error('File read unsuccessful')
    end
    
    img(:,:,j) = reshape(image,sizexy);
    j = j + 1;
end

%     imagesc(image)
%     head = [head,'file: ',fname,' image #',num2str(i)];
end