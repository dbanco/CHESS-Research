function [polar_image, rad_domain, theta_domain,theta_ind] = extract_ring( image, r1, r2,center,dtheta,num_theta)
%extract_ring Summary of this function goes here
%   Detailed explanation goes here

rad_domain = r1:r2;
num_rad = numel(rad_domain);
polar_image = zeros(num_rad,num_theta);
theta_domain  = linspace(0,2*pi-dtheta,num_theta);

for i = 1:num_rad
    [projection, theta_domain,theta_ind] = azimuthal_projection(image,center,rad_domain(i),theta_domain);
    polar_image(i,theta_ind) = projection;
end

rowSum = sum(polar_image,1);

theta_ind = logical(rowSum > 0);
theta_domain(:,rowSum == 0) = [];
polar_image(:,rowSum == 0) = [];


