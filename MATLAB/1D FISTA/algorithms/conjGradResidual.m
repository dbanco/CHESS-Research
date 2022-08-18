function R = conjGradResidual(A0ft_stack,X,Xft,AtB)
%conjGradResidual Compute residual for conjugate gradient that includes 
% difference matrix
%   Detailed explanation goes here

T = size(B,2);
[N,M] = size(A0ft_stack);
Y = zeros(N,M,T);
Diff1 = cell(N,M,T-1);
Y2 = cell(N,M,T);

for t = 1:(T-1)
    Diff1(:,:,t) = X{t+1}-X{t};
end

Y2(:,:,1) = -Diff1(:,:,1);
Y2(:,:,T) =  Diff1(:,:,T);
for t= 2:(T-1)
    Y2(:,:,t) = Diff1(:,:,t-1) - Diff1(:,:,t);
end

for t = 1:T
    Y(:,:,t) = AtR_ft_1D(A0ft_stack,Ax_ft_1D(A0ft_stack,Xft{t}));
end

R = Y + Y2 + AtB;
