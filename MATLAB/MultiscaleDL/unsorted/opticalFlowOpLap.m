function out = opticalFlowOpLap(x,K,trans)
%opticalFlowOp

if nargin < 5
    trans = 0;
end

x = squeeze(x);
U = size(x,2)/K;

Fx = zeros(size(x));
Ft = zeros(size(x));

% Apply space time differences to each of K atoms
for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;

    Fx(:,i1:i2,:) = laplaceFilt(x(:,i1:i2,:));
    Ft(:,i1:i2,:) = timeDiff(x(:,i1:i2,:));
end

% Option to additionally apply operation transposed
if trans
    Ax = Fx + Ft;
    for k = 1:K
        i1 = 1+(k-1)*U;
        i2 = k*U;
    
        Fx(:,i1:i2,:) = laplaceFilt(Ax(:,i1:i2,:),1);
        Ft(:,i1:i2,:) = timeDiff(Ax(:,i1:i2,:),1);
    end
    out = Fx + Ft;
else
    out = Fx + Ft;
end

end