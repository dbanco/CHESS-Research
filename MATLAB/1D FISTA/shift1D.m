function [ image ] = shift1D( image, rows)
%shift2D circularly shifts image by row, col pixels
% Inputs:
% image - a vector
% rows - any integer number of rows
% Output:
% image - circularly shifted image
try
    for i = 1:rows
       image = [image(end) ; image(1:end-1)]; 
    end    
catch
    for i = 1:rows
       image = [image(end) , image(1:end-1)]; 
    end 
end
end
