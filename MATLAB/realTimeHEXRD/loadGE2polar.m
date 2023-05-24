function [b, rad, az] = loadGE2polar(Trange,fname,center,r1,r2)

T = numel(Trange);
num_rad = r2-r1+1;
num_eta = round(2*pi*r2);
b = zeros(num_rad,num_eta,T);
rad = zeros(num_rad,T);
az = zeros(num_eta,T);

img = readGE2img(fname,Trange);

dtheta = 2*pi./num_eta;
for t = 1:T
    [ring, radt, azt] = extract_ring( img(:,:,t),r1,r2,center,dtheta,num_eta );
    b(:,:,t) = ring;
    rad(:,t) = radt;
    az(:,t) = azt;
end



end