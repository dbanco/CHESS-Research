function A0_stack = unshifted_basis_matrix_stack(P)
%unshifted_basis_matrix_stack Generates many zero mean gaussian  
% basis function matrices that sum to 1 using provided parameters
%
% Inputs:
% P:
% var_theta -vector of theta variances
% var_rad - vector of radial variances
% dtheta - difference in theta between each pixel
% drad - difference in radius between each pixel
% num_theta - image size in theta direction
% num_rad - image size in radial direction
% weight - 1 or 0, to weight or not to weight
% alphap - total penalty
% betap  - variance based penalty
%
% Outputs:
% A0_stack - (n x m x t x r) array
%             n = num_theta
%             m = num_rad
%             t = numel(var_theta)
%             r = numel(var_rad)

A0_stack = zeros(P.num_rad,P.num_theta,numel(P.var_theta),numel(P.var_rad));
for t = 1:numel(P.var_theta)
    for r = 1:numel(P.var_rad)
        A0 = gaussian_basis_wrap_2D(P.num_theta,P.dtheta,  0,  P.var_theta(t),...
                                    P.num_rad,  P.drad,    0,  P.var_rad(r));
        A0_stack(:,:,t,r) = P.alphap*2.^(-P.weight*P.betap*...
                                          P.dtheta^2/P.var_theta(t))*...
                                          P.drad^2/P.var_rad(r)*A0;
    end
end
end

