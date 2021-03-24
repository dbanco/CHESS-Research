function B = gaussian_basis_wrap_1D_norm2(num_theta,mean_theta,var_theta)
%gaussian_basis_wrap_2D Generates gaussian basis function matrix
% normalized to sum to 1
%
% Inputs:
% num_theta - image size in theta direction
% dtheta - difference in theta between each pixel
% mean_theta - mean of gaussian basis function in theta
% var_theta - variance of gaussian basis function in theta
% (n) is the size of the image and basis functions

% Outputs:
% Ax - (n x m) array

% Compute theta distances with wrapping
idx = 1:num_theta;
wrapN = @(x, N) (1 + mod(x-1, N));
opposite = (idx(wrapN(floor(mean_theta-num_theta/2),num_theta)) +... 
            idx(wrapN(ceil(mean_theta-num_theta/2),num_theta)))/2;
if opposite == mean_theta
    opposite = 0.5;
end
dist1 = abs(mean_theta - idx);
dist2 = num_theta/2 - abs(opposite - idx);
dist = min(dist1,dist2);
dist_sq_theta = dist.^2;    % num_theta length vector

% Compute values
B = exp(-dist_sq_theta/(2*var_theta) )';
B = B/norm(B(:));
end

