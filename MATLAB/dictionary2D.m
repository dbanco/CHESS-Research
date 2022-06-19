function D = dictionary2D(P)
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

K = size(P.stds,1);
D = zeros(P.Ny,P.Nx,K);
try
    Pmeans = P.means;
catch
    Pmeans = ones(K,2);
end
for i = 1:K
    d = gaussian_basis_wrap_2D(P.Ny, Pmeans(i,1), P.stds(i,1),...
                               P.Nx, Pmeans(i,2), P.stds(i,2),'2-norm');               
    D(:,:,i) = d;
end

end

