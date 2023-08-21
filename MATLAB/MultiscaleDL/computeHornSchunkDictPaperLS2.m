function [u,v,Fx,Fy,Ft] = computeHornSchunkDictPaperLS2(data,K,Uvel,Vvel,smoothness,maxIters)
%computeHornSchunkDict


if nargin <5
    smoothness = 1;
end
if nargin <6
    maxIters = 10;
end

data = squeeze(data);
KU = size(data,2);
U = KU/K;

u = zeros(size(data));
v = zeros(size(data));
Fx = zeros(size(data));
Fy = zeros(size(data));
Ft = zeros(size(data));

for k = 1:K
    i1 = 1+(k-1)*U;
    i2 = k*U;
    if isempty(Uvel)
        u0 = zeros(size(data(:,i1:i2,:)));
        v0 = u0;
    else
        u0 = Uvel(:,i1:i2,:);
        v0 = Vvel(:,i1:i2,:);
    end
    [uk,vk,Fxk,Fyk,Ftk] = computeHornSchunkPaperLS2(data(:,i1:i2,:),u0,v0,smoothness,maxIters);
    
    u(:,i1:i2,:) = uk;
    v(:,i1:i2,:) = vk;
    Fx(:,i1:i2,:) = Fxk;
    Fy(:,i1:i2,:) = Fyk;
    Ft(:,i1:i2,:) = Ftk;
end

end