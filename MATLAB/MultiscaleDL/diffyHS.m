function d = diffyHS(x)
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

kernel = zeros(2,2,2);
kernel(:,:,1) = [ -1 -1; 1 1];
kernel(:,:,2) = [ -1 -1; 1 1];

% d = convn(x,rot90(kernel,2),'valid');
d = convn(x,kernel,'valid');
end
