function A0_stack = unshifted_basis_matrix_ft_stack_norm2_3D(P)
%unshifted_basis_matrix_stack_norm Generates fft2 of many zero mean gaussian  
% basis function matrices that sum to 1 using provided parameters
%
% Inputs:
% P:
% var_theta -vector of theta variances
% var_rad - vector of radial variances
% var_omega-vector of omega variances
% dtheta - difference in theta between each pixel
% drad - difference in radius between each pixel
% domega - difference in omega between each pixel
% num_theta - image size in theta direction
% num_rad - image size in radial direction
% num_omega - image size in omega direction
% weight - 1 or 0, to weight or not to weight
% 
%
% Outputs:
% A0ft_stack - (n x m x p x t x r x w) array
%             n = num_theta
%             m = num_rad
%             p = num_omega
%             t = numel(var_theta)
%             r = numel(var_rad)
%             w = numel(var_omega)
A0_stack = zeros(P.num_rad,P.num_theta,P.num_omega,numel(P.var_theta),numel(P.var_rad),numel(P.var_omega));
for t = 1:numel(P.var_theta)
    for r = 1:numel(P.var_rad)
        for w = 1:numel(P.var_omega)
            A0 = gaussian_basis_wrap_3D_norm2(P.num_theta,P.dtheta, 0, P.var_theta(t),...
                                              P.num_rad,  P.drad,   0, P.var_rad(r),...
                                              P.num_omega,P.domega, 0, P.var_omega(w));
            A0_stack(:,:,:,t,r,w) = fftn(A0);
        end
    end
end
end

