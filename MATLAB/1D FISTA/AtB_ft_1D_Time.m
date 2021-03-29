function AtB = AtB_ft_1D_Time(A0ft_stack,B,Bnorms)
%conjGradResidual Compute residual for conjugate gradient that includes 
% difference matrix
%   Detailed explanation goes here

[N,M] = size(A0ft_stack);
T = size(B,2);
AtB = zeros(N,M,T);
for t = 1:T
    AtB(:,:,t) = AtR_ft_1D(A0ft_stack/Bnorms(t),B(:,t));
end