function A0_stack = unshifted_basis_matrix_stack(P)
%unshifted_basis_matrix_stack Generates many zero mean gaussian  
% basis function matrices that sum to 1 using provided parameters
%
% Inputs:
% P:
% var_theta -vector of theta variances
% var_rad - vector of radial variances
% num_theta - image size in theta direction
% num_rad - image size in radial direction

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
        switch P.basis
            case 'norm2'
                A0 = gaussian_basis_wrap_2D_norm2(P.num_theta,1,  0, P.var_theta(t),...
                                                  P.num_rad,  1,  0, P.var_rad(r));
            case 'max'
                A0 = gaussian_basis_wrap_2D(P.num_theta, 0, P.var_theta(t),...
                                            P.num_rad,   0, P.var_rad(r),'max')';
            case 'norm1'
                A0 = gaussian_basis_wrap_2D_norm1(P.num_theta,1,  0, P.var_theta(t),...
                                                  P.num_rad,  1,  0, P.var_rad(r));
            case 'norm'
                A0 = gaussian_basis_wrap_2D_norm(P.num_theta,1,  0, P.var_theta(t),...
                                                  P.num_rad,  1,  0, P.var_rad(r));
        end
        A0_stack(:,:,t,r) = A0;
    end
end
end
