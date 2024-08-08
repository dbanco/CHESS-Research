function [Ft, addTerm] = timeDiffDict(X,K,U,trans)
%timeDiffDict 

if nargin < 4
    trans = 0;
end

X = squeeze(X);
Ft = zeros(size(X));
for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;
    Ft(:,i1:i2,:) = -timeDiff(X(:,i1:i2,:),trans);
end

if nargout == 2
    addTerm = X;
    addTerm(Ft==0) = 0;
end

end