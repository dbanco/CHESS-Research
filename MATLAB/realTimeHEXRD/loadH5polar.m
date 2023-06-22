function [b, rad, az] = loadH5polar(t1,t2,fname,center,r1,r2)

T = t2-t1+1;
num_rad = r2-r1+1;
num_eta = round(2*pi*r2);
dtheta = 2*pi./num_eta;
b = zeros(num_rad,num_eta,T);
rad = zeros(num_rad,T);
az = zeros(num_eta,T);
img = squeeze(h5read( fname,'/imageseries/images', ...
                      [1 1 t1],[3072 3888 t2] ));
for t = 1:T
    [ring, radt, azt,ind] = extract_ring( img(:,:,t),r1,r2,center,dtheta,num_eta );
    b(:,ind,t) = ring;
    rad(:,t) = radt;
    az(ind,t) = azt;
end
rowSum = sum(b,[1,3]);
b(:,rowSum == 0,:) = [];
az(rowSum == 0,:) = [];
end