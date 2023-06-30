function [uFx,vFy,Ft] = diffDictPaperTrans(X,K,U,u,v)
%sobelDictTrans

X = squeeze(X);
vFy = zeros(size(X));
uFx = zeros(size(X));
for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;
    dataPad = padarray(X(:,i1:i2,:).*v(:,i1:i2,:),[1 1 1],0,'pre');
    vFy(:,i1:i2,:) = diffyPaper(dataPad);
    dataPad = padarray(X(:,i1:i2,:).*u(:,i1:i2,:),[1 1 1],0,'pre');
    uFx(:,i1:i2,:) = diffxPaper(dataPad);
    dataPad = padarray(X(:,i1:i2,:),[1 1 1],0,'pre');
    Ft(:,i1:i2,:) = -difftPaper(dataPad);
end

end