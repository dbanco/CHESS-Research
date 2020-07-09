function [ image ] = shift1D_array( image, rows)
%shift2D circularly shifts image by row, col pixels
% Inputs:
% image - an array
% rows - any integer number of rows
% Output:
% image - circularly shifted image
N = size(image,1);
image = [image((N-rows+1):N,:) ; image(1:N-rows,:)]; 

end
