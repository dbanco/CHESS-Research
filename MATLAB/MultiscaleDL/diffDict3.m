function [Fx,Fy,Ft] = diffDict3(X,K,U)
%sobelDict 

X = squeeze(X);
Fy = zeros(size(X));
Fx = zeros(size(X));

for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;
    dataPad = padarray(X(:,i1:i2,:),[1 1],0,'both');
    Fy(:,i1:i2,:) = diffyHS3(dataPad);
    Fx(:,i1:i2,:) = diffxHS3(dataPad);
    dataPad = padarray(data,[1 1 1],0,'pre');
    Ft(:,i1:i2,:) = difftHS2(dataPad);
end

end