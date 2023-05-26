function [b, rad, az] = loadGE2polar(Trange,fname,center,r1,r2)

T = numel(Trange);
num_rad = r2-r1+1;
num_eta = round(2*pi*r2);
b = zeros(num_rad,num_eta,T);
rad = zeros(num_rad,T);
az = zeros(num_eta,T);

img = readGE2img(fname,Trange);
img = img(1:1000,:);
dtheta = 2*pi./num_eta;
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