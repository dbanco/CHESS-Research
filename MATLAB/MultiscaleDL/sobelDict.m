function [Fx,Fy] = sobelDict(X,K,U,trans)
%sobelDict 

if nargin < 4
    trans = 0;
end
if trans

else
X = squeeze(X);
Fy = zeros(size(X));
Fx = zeros(size(X));
for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;
    Fy(:,i1:i2,:) = sobel(X(:,i1:i2,:),1,trans);
    Fx(:,i1:i2,:) = sobel(X(:,i1:i2,:),2,trans);
end

end