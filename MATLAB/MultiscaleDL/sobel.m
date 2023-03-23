function d = sobel(x,dim,trans)
%sobel Sobel edge detector in row dimenion(1) or column dimension(2)
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

kernel = [-1 -2 -1; 0 0 0; 1 2 1];

if nargin < 3
    trans = 0;
end

if dim == 2
    kernel = kernel';
elseif dim ~= 1
    error('Invalid dimension')
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
