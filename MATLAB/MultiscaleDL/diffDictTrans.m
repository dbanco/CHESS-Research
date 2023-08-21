function [Fx,Fy,Ft] = diffDictTrans(X,K,U,u,v)
%sobelDict 

X = squeeze(X);
Fy = zeros(size(X));
Fx = zeros(size(X));

for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;
    dataPad = padarray(-X(:,i1:i2,:).*u(:,i1:i2,:),[1 1 1],0,'post');
    Fx(:,i1:i2,:) = diffxHS(dataPad);
    dataPad = padarray(-X(:,i1:i2,:).*v(:,i1:i2,:),[1 1 1],0,'post');
    Fy(:,i1:i2,:) = diffyHS(dataPad);
    dataPad = padarray(-X(:,i1:i2,:),[1 1 1],0,'post');
    Ft(:,i1:i2,:) = difftHS(dataPad);
end

end