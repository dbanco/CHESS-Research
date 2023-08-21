function [Fx,Fy,Ft] = diffDict(X,K,U)
%sobelDict 

X = squeeze(X);
Fy = zeros(size(X));
Fx = zeros(size(X));

for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;
    dataPad = padarray(X(:,i1:i2,:),[1 1 1],0,'pre');
    Fy(:,i1:i2,:) = diffyHS(dataPad);
    Fx(:,i1:i2,:) = diffxHS(dataPad);
    Ft(:,i1:i2,:) = difftHS(dataPad);
end

end