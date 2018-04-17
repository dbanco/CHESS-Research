function B = gaussian_basis_wrap_2D_norm(num_theta,dtheta,mean_theta,var_theta,...
                                         num_rad,drad,mean_rad,var_rad )
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

% Compute theta distances with wrapping
idx = 1:num_theta;
wrapN = @(x, N) (1 + mod(x-1, N));
opposite = (idx(wrapN(mean_theta-floor(num_theta/2),num_theta)+1) +... 
            idx(wrapN(mean_theta-ceil(num_theta/2),num_theta)+1))/2;
if opposite == mean_theta
    opposite = 0.5;
end
dist1 = abs(mean_theta - idx);
dist2 = num_theta/2 - abs(opposite - idx);
dist = min(dist1,dist2).*dtheta;
dist_sq_theta = dist.^2;    % num_theta length vector

% Compute radial distance with wrapping
idx = 1:num_rad;
opposite = (idx(wrapN(mean_rad-floor(num_rad/2),num_rad)+1) +...
            idx(wrapN(mean_rad-ceil(num_rad/2),num_rad)+1))/2;
if opposite == mean_rad
    opposite = 0.5;
end
dist1 = abs(mean_rad - idx);
dist2 = num_rad/2 - abs(opposite - idx);
dist = min(dist1,dist2).*drad;
dist_sq_rad = dist.^2;    % num_r length vector

% Create matrix of radial sqaured distances and matrix of theta squared 
% distances to compute gaussian basis function value at each point
all_dist_sq_theta = zeros(num_rad,num_theta);
all_dist_sq_rad = zeros(num_rad,num_theta);
for i = 1:num_rad
    all_dist_sq_theta(i,:) = dist_sq_theta';
    all_dist_sq_rad(i,:) = dist_sq_rad(i).*ones(num_theta,1);
end

% Compute values
B = exp(-all_dist_sq_theta/(2*var_theta) - all_dist_sq_rad/(2*var_rad));
B = B/norm(B(:));
end

