function [ image ] = shift2D( image, rows, cols )
%shift2D shifts image by row, col pixels

for i = 1:rows
   image = [image(end,:) ; image(1:end-1,:)]; 
end    
    
for j = 1:cols
   image = [image(:,end) , image(:,1:end-1)]; 
end

end
