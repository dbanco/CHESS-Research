function removeWhiteSpace(filePath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Read the image file
img = imread(filePath);

% Convert to grayscale and invert colors
gray = 255 * (rgb2gray(img) < 128);

% Find non-zero pixels (image content)
[row,col] = find(gray);

% Crop the image using the bounding box
cropped = img( min(row):max(row), min(col):max(col),:);

% Save the cropped image overwriting
imwrite(cropped, filePath);
end