function [dd,center]= dSampleCenter2(d, m, center)
% Downsample with stride m, centered around `center`
% using symmetric padding of length m

if m == 1
    dd = d;
else
    d = padarray(d,[0,m-1,0],'both');
    center = center + m - 1;
    N = size(d,2);
    inds = [fliplr(center-m:-m:1),center:m:N];
    d1 = lowpassM(d,m,0);
    dd = d1(inds);
    center = find(inds==center);
end

end
