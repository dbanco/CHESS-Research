function B = gaussian_basis_2D(num_theta,  mean_theta,var_theta, num_rad,  mean_rad,var_rad)
%gaussian_basis_wrap_2D Generates gaussian basis function matrix
% normalized to sum to 1
%
% Inputs:
% num_theta - image size in theta direction
% dtheta - difference in theta between each pixel
% mean_theta - mean of gaussian basis function in theta
% var_theta - variance of gaussian basis function in theta
% num_rad - image size in radial direction
% drad - difference in radius between each pixel
% mean_rad - radial mean of gaussian basis function
% var_rad - radial variance of gaussian basis function
% (n x m) is the size of the image and basis functions
% (t x r) indexes the basis function by theta variance and radial variance
%
% Outputs:
% Ax - (n x m) array

% Compute theta distances without wrapping
[ind_r,ind_t] = find(ones(num_rad,num_theta));
dist_sq_theta = (ind_t - mean_theta).^2; 
dist_sq_rad = (ind_r - mean_rad).^2;

% Compute values
B = exp(-dist_sq_theta/(2*var_theta) - dist_sq_rad/(2*var_rad));
B = B/sum(B(:));
B = reshape(B,[num_rad,num_theta]);
end

