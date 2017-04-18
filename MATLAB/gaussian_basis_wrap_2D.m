function B = gaussian_basis_wrap_2D( num_theta,dtheta,mean_theta,var_theta,...
                                     num_rad,drad,mean_rad,var_rad )
%gaussian_basis_wrap_2D Generates basis function matrix
%   Detailed explanation goes here

% Compute theta distances with wrapping
idx = 1:num_theta;
wrapN = @(x, N) (1 + mod(x-1, N));
opposite = (idx(wrapN(mean_theta-floor(num_theta/2),num_theta)+1) +... 
            idx(wrapN(mean_theta-ceil(num_theta/2),num_theta)+1))/2;
dist1 = abs(mean_theta - idx);
dist2 = num_theta/2 - abs(opposite - idx);
dist = min(dist1,dist2).*dtheta;
dist_sq_theta = dist.^2;    % num_theta length vector

% Compute radial distance with wrapping
idx = 1:num_rad;
opposite = (idx(wrapN(mean_rad-floor(num_rad/2),num_rad)+1) +...
            idx(wrapN(mean_rad-ceil(num_rad/2),num_rad)+1))/2;
dist1 = abs(mean_rad - idx);
dist2 = num_rad/2 - abs(opposite - idx);
dist = min(dist1,dist2).*drad;
dist_sq_rad = dist.^2;    % num_r length vector

all_dist_sq_theta = dist_sq_theta';
all_dist_sq_rad = dist_sq_rad(1).*ones(num_theta,1);
for i = 2:num_rad
    all_dist_sq_theta = [all_dist_sq_theta; dist_sq_theta'];
    all_dist_sq_rad = [all_dist_sq_rad; dist_sq_rad(i).*ones(num_theta,1)];
end
% Build basis vector
B = exp(-all_dist_sq_theta./var_theta - all_dist_sq_rad./var_rad);
B = reshape(B,num_theta,num_rad)';
end

