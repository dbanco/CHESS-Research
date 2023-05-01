function d = diffy(x,trans)
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

if nargin < 2
    trans = 0;
end

if trans == 0
    d = x;
    d(2:end,:,:) = x(2:end,:,:)-x(1:end-1,:,:);
elseif trans == 1
%     xpad = padarray(x,[1,1],0);
%     k = zeros(size(xpad,1),size(xpad,2));
%     k(1:3,1:3) = kernel;
%     k = circshift(k,[-1,-1]);
%     xft = fft2(xpad);
%     kft = fft2(k);
%     d = ifft2(bsxfun(@times,conj(kft),xft),'symmetric');
%     d = d(2:end-1,2:end-1,:);
end

end
