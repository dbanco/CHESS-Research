function A0_stack = unshifted_basis_vector_ft_stack_zpad(P)
%unshifted_basis_matrix_stack_norm Generates fft2 of many zero mean gaussian  
% basis function matrices that sum to 1 using provided parameters
%
% Inputs:
% P:
% var_theta -vector of theta variances
% num_theta - image size in theta direction
%
% zPad
%
% Outputs:
% A0ft_stack - (n x t) array
%             n = num_theta
%             t = numel(var_theta)
zPad = P.params.zeroPad;
if sum(zPad)
    A0_stack = zeros(P.num_theta + 2*zPad, numel(P.var_theta));
    for t = 1:numel(P.var_theta)
        switch P.basis
            case 'norm2'
                A0 = gaussian_basis_wrap_1D_norm2(P.num_theta, 1, P.var_theta(t));
        end
        A0_stack(:,t) = fft(zeroPad(A0,zPad));
    end
else
    A0_stack = zeros(P.num_theta,numel(P.var_theta));
    for t = 1:numel(P.var_theta)
        switch P.basis
            case 'norm2'
                A0 = gaussian_basis_wrap_1D_norm2(P.num_theta, 1, P.var_theta(t));
        end                     
        A0_stack(:,t) = fft(A0);
    end
end
end

