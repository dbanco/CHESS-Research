function [uFx,vFy] = sobelDictTransFix(X,K,U,u,v)
%sobelDictTrans

X = squeeze(X);
vFy = zeros(size(X));
uFx = zeros(size(X));
for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;
    vFy(:,i1:i2,:) = sobel(X(:,i1:i2,:).*v(:,i1:i2,:),1);
    uFx(:,i1:i2,:) = sobel(X(:,i1:i2,:).*u(:,i1:i2,:),2);
end

end