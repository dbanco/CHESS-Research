function [du, center] = uSampleCenter2(d, m, center)
% Upsample signal `d` (1D row) by stride m, centered around `center`
% into output of length N, using symmetric padding of m on both sides

if m == 1
    du = d;
else
    N = size(d,2);
    du = zeros(1,N*m);
    inds = [fliplr(center-m:-m:1),center:m:N*m];
    du(1,inds) = d;
    preZeros = (m-1) - (inds(1)-1); % # zeros needed - # zeros we have
    postZeros = (m-1) - (N*m-inds(end));% # zeros needed - # zeros we have
    du = padarray(du,[0,preZeros],0,'pre');
    du = padarray(du,[0,postZeros],0,'post');
    center = center + preZeros;
    du = lowpassM(du,m,0);
    
end


end