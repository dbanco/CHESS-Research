function Y = conjGradAxDx(A0ft_stack,X,params)
%conjGradResidual Compute residual for conjugate gradient that includes 
% difference matrix
%   Detailed explanation goes here

T = size(X,3);
[N,M] = size(A0ft_stack);

Y1 = zeros(N,M,T);
for t = 1:T
    Y1(:,:,t) = AtR_ft_1D(A0ft_stack,Ax_ft_1D(A0ft_stack,X(:,:,t)));
end

Y2 = zeros(N,M,T);
Diff1 = zeros(N,M,T-1);
for t = 1:(T-1)
    Diff1(:,:,t) = X(:,:,t+1)-X(:,:,t);
end
Y2(:,:,1) = -Diff1(:,:,1);
Y2(:,:,T) =  Diff1(:,:,T-1);
for t= 2:(T-1)
    Y2(:,:,t) = Diff1(:,:,t-1) - Diff1(:,:,t);
end

Y = Y1 + params.rho*Y2;
