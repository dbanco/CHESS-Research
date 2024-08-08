function Aty = opticalFlowOpTran(y,u,v,K)
%opticalFlowOpTran
[N1,N2,KJ,T] = size(y);
y = squeeze(y);
U = size(y,2)/K;

Fx = zeros(size(y));
Fy = zeros(size(y));
Ft = zeros(size(y));
for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;

    Fy(:,i1:i2,:) = sobel(y(:,i1:i2,:),1,1);
    Fx(:,i1:i2,:) = sobel(y(:,i1:i2,:),2,1);
    Ft(:,i1:i2,:) = timeDiff(y(:,i1:i2,:),1);
end
Aty = Fx.*u + Fy.*v + Ft;
Aty = reshape(Aty,[N1,N2,KJ,T]);
end