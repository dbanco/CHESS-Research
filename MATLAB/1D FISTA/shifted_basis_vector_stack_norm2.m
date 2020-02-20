function A_stack = shifted_basis_vector_stack_norm2(P)
%unshifted_basis_matrix_stack Generates many zero mean gaussian  
% basis function matrices that sum to 1 using provided parameters
%
% Inputs:
% P:
% var_theta -vector of theta variances
% dtheta - difference in theta between each pixel
% num_theta - image size in theta direction
%
% Outputs:
% A0ft_stack - (n x t) array
%             n = num_theta
%             t = numel(var_theta)
A_stack = zeros(P.num_theta,numel(P.var_theta)*P.num_theta);
N = P.num_theta;
for t = 1:numel(P.var_theta)
    A0 = gaussian_basis_wrap_1D_norm2(P.num_theta, 0, P.var_theta(t));                      
    for i = 1:N
        A_stack(:,(t-1)*N + i) = shift1D(A0,i);
    end
end


end
