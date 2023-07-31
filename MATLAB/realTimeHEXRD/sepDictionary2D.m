function [sepD1,sepD2] = sepDictionary2D(P)
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


sepD1 = zeros(P.N1,1,P.K1,1);
sepD2 = zeros(1,P.N2,1,P.K2);
% try
%     Pmeans = P.means;
% catch
%     Pmeans = ones(P.K,2);
% end
for i1 = 1:P.K1
    for i2 = 1:P.K2
        d = gaussian_basis_wrap_2D(P.N1, P.mu1, P.sigma1(i1),...
                                   P.N2, P.mu2, P.sigma2(i2),'2-norm');    
        [U,~,V] = svd(d);
        
        sepD1(:,1,i1) = abs(U(:,1));
        sepD2 (1,:,i2)= abs(V(:,1));
    end
end

end

