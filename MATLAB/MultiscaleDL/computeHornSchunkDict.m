function [u,v,Fy,Fx,Ft] = computeHornSchunkDict(data,K,smoothness,maxIters)
%computeHornSchunkDict

if nargin <3
    smoothness = 1;
end
if nargin <4
    maxIters = 10;
end

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
    [uk,vk,Fyk,Fxk,Ftk] = computeHornSchunk(data(:,i1:i2,:),smoothness,maxIters);
    u(:,i1:i2,:) = uk;
    v(:,i1:i2,:) = vk;
    Fx(:,i1:i2,:) = Fxk;
    Fy(:,i1:i2,:) = Fyk;
    Ft(:,i1:i2,:) = Ftk;
end

norm(Fx(:))
norm(Fy(:))
norm(Ft(:))
norm(u(:))
norm(v(:))

end