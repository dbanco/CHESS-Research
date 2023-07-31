function [b, rad, az] = loadH5polarDex(t1,t2,fname,fname2,center,r1,r2)

T = t2-t1+1;
num_rad = r2-r1+1;
num_eta = round(2*pi*r2);
dtheta = 2*pi./num_eta;
b = zeros(num_rad,num_eta,1);
rad = zeros(num_rad,1);
az = zeros(num_eta,1);
img1 = squeeze(h5read( fname,'/imageseries/images', ...
                      [1 1 t1],[3072 3888 T] ));
img2 = squeeze(h5read( fname2,'/imageseries/images', ...
                      [1 1 t1],[3072 3888 T] ));
img1 = squeeze(max(img1,[],3));
img2 = squeeze(max(img2,[],3));
gap = round(10.5/0.0748);
full_dex = [flip(img1',1),zeros(3888,gap),flip(img2',2)];
                  
for t = 1:size(full_dex,3)
    [ring, radt, azt,ind] = extract_ring( full_dex(:,:,t),r1,r2,center,dtheta,num_eta );
    b(:,ind,t) = ring;
    rad(:,t) = radt;
    az(ind,t) = azt;
end
rowSum = sum(b,[1,3]);
b(:,rowSum == 0,:) = [];
az(rowSum == 0,:) = [];
end