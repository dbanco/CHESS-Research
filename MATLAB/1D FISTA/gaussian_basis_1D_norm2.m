function B = gaussian_basis_1D_norm2(num_theta,  mean_theta,var_theta)
%gaussian_basis_wrap_2D Generates gaussian basis function matrix
% normalized to sum to 1
%
% Inputs:
% num_theta - image size in theta direction
% dtheta - difference in theta between each pixel
% mean_theta - mean of gaussian basis function in theta
% var_theta - variance of gaussian basis function in theta
% num_rad - image size in radial direction
% (n x m) is the size of the image and basis functions
% (t x r) indexes the basis function by theta variance and radial variance
%
% Outputs:
% Ax - (n x m) array

% Compute theta distances without wrapping
[ind_t] = find(ones(num_theta,1));
dist_sq_theta = (ind_t - mean_theta).^2; 

% Compute values
B = exp(-dist_sq_theta/(2*var_theta));
B = B/norm(B(:));
end

