function A0_stack = unshifted_basis_matrix_stack_norm(var_theta,var_rad,dtheta,drad,num_theta,num_rad)
%unshifted_basis_matrix_stack Summary of this function goes here
%   Detailed explanation goes here

A0_stack = zeros(num_rad,num_theta,numel(var_theta),numel(var_rad));
for t = 1:numel(var_theta)
    for r = 1:numel(var_rad)
        A0 = gaussian_basis_wrap_2D_norm(num_theta,dtheta,  0,  var_theta(t),...
                                         num_rad,  drad,    0,  var_rad(r));
        A0_stack(:,:,t,r) = A0;
    end
end
end

