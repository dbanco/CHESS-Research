function Y = Ax_ft_1D_Time(A0ft_stack,X)
%conjGradResidual Compute residual for conjugate gradient that includes 
% difference matrix
%   Detailed explanation goes here

T = size(X,3);
[N,M] = size(A0ft_stack);

Y = zeros(N,T);
for t = 1:T
    Y(:,t) = Ax_ft_1D(A0ft_stack,X(:,:,t));
end
