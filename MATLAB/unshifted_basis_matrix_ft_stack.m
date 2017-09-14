function A0ft_stack = unshifted_basis_matrix_ft_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad)
%unshifted_basis_matrix_ft_stack Generates fft2 of many zero mean gaussian  
% basis function matrices with a max value of 1 using provided parameters
%
% Inputs:
% var_theta -vector of theta variances
% var_rad - vector of radial variances
% dtheta - difference in theta between each pixel
% drad - difference in radius between each pixel
% num_theta - image size in theta direction
% num_rad - image size in radial direction
%
% Outputs:
% A0ft_stack - (n x m x t x r) array
%             n = num_theta
%             m = num_rad
%             t = numel(var_theta)
%             r = numel(var_rad)

A0ft_stack = zeros(num_rad,num_theta,numel(var_theta),numel(var_rad));
for t = 1:numel(var_theta)
    for r = 1:numel(var_rad)
        A0 = gaussian_basis_wrap_2D(num_theta,dtheta,  0,  var_theta(t),...
                                    num_rad,  drad,    0,  var_rad(r));
        A0ft_stack(:,:,t,r) = fft2(A0);
    end
end
end
