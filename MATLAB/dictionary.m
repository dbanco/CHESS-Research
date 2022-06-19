function A0_stack = dictionary(P)
%dictionary Generates zero mean gaussian  
% basis function vectors with unit 2-norm
%
% Inputs:
% P.var_theta -vector of theta variances
% P.dtheta - difference in theta between each pixel
% P.num_theta - image size in theta direction 
%
% Outputs:
% A0ft_stack - Dictionary atoms [N,K]

K = numel(P.stds);
A0_stack = zeros(P.N,K);
try
    Pmeans = P.means;
catch
    Pmeans = ones(K,1);
end
for i = 1:K
    A0 = gaussian_basis_wrap_1D(P.N, Pmeans(i), P.stds(i),'2-norm');               
    A0_stack(:,i) = A0;
end

end

