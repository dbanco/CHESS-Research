function [polar_image, rad_domain, theta_domain] = extract_ring( image, r1, r2,center,dtheta,num_theta)
%extract_ring Summary of this function goes here
%   Detailed explanation goes here

rad_domain = r1:r2;
num_rad = numel(rad_domain);
polar_image = zeros(num_rad,num_theta);

for i = 1:num_rad
    [projection, theta_domain] = azimuthal_projection(image,center,rad_domain(i),0,2*pi-dtheta,num_theta);
    polar_image(i,:) = projection;
end

