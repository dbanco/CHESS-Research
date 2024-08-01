function d = laplaceFilt(x,trans)
%laplace Laplace smoothing filt in row dimenion(1) or column dimension(2)
%Equivalent to zero padding x and multiplying by BCCB of kernel
% xpad = padarray(x,[1,1],0);
% k = zeros(size(xpad,1),size(xpad,2));
% k(1:3,1:3) = kernel;
% k = circshift(k,[-1,-1]);
% imagesc(k)
% xft = fft2(xpad);
% kft = fft2(k);
% dx = ifft2(xft.*kft,'symmetric');
% dx = dx2(2:end-1,2:end-1,:);

kernel = [-1 -1 -1; -1 8 -1; -1 -1 -1]/8;

if nargin < 2
    trans = 0;
end

if trans == 0
    d = convn(x,kernel,'same');
elseif trans == 1
    xpad = padarray(x,[1,1],0);
    k = zeros(size(xpad,1),size(xpad,2));
    k(1:3,1:3) = kernel;
    k = circshift(k,[-1,-1]);
    xft = fft2(xpad);
    kft = fft2(k);
    d = ifft2(bsxfun(@times,conj(kft),xft),'symmetric');
    d = d(2:end-1,2:end-1,:);
end

end
