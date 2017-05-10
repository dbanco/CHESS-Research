function [ image ] = shift2D( image, rows, cols )
%shift2D circularly shifts image by row, col pixels
% Inputs:
% image - any size image
% rows - any integer number of rows
% cols - any integer number of cols
% Output:
% image - circularly shifted image

for i = 1:rows
   image = [image(end,:) ; image(1:end-1,:)]; 
end    
    
for j = 1:cols
   image = [image(:,end) , image(:,1:end-1)]; 
end

end
