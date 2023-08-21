function [u,v,Fx,Fy,Ft] = computeHornSchunkPaperLS3(data,u0,v0,smoothness,maxIters,tol)
%computeHornSchunk Summary of this function goes here
%   Detailed explanation goes here

if nargin <4 || isempty(smoothness)
    smoothness = 1;
end
if nargin <5 || isempty(maxIters)
    maxIters = 100;
end
if nargin < 6 || isempty(tol)
  tol = 1e-6;
end

% Compute partial derivatives
dataPad = gpuArray(padarray(data,[1 1],0,'both'));
Fx = diffxHS3(dataPad);
Fy = diffyHS3(dataPad);
dataPad = gpuArray(padarray(data,[1 1 1],0,'pre'));
Ft = difftHS2(dataPad);

% Set up linear system
b = [-Fx(:).*Ft(:); -Fy(:).*Ft(:)];
N = numel(data);
Aop = @(uv) Aoperator(Fx(:),Fy(:),uv,N,smoothness,size(u0));

% Solve with conjugate gradient
uv0 = [u0(:);v0(:)];
[uv,~] = pcg(@(u) Aop(u), b(:),tol,maxIters,[],[],uv0);
u = reshape(uv(1:N),    size(u0));
v = reshape(uv(N+1:end),size(v0));
end

function Ax = Aoperator(Fx,Fy,uv,N,smoothness,uSize)
    u = uv(1:N);
    v = uv((N+1):end);
    FxFy = Fx.*Fy;
    [nablaU2,nablaV2] = laplaceOp2(reshape(u,uSize),reshape(v,uSize));
%     Au = smoothness*2*(diffx(diffx(u),1) + diffy(diffy(u),1)) + Bu;
%     Av = smoothness*2*(diffx(diffx(v),1) + diffy(diffy(v),1)) + Bv;
    Au = Fx.^2.*u + FxFy.*v - smoothness*nablaU2(:) ;
    Av = Fy.^2.*v + FxFy.*u - smoothness*nablaV2(:) ;
    Ax = [Au; Av];
end


% function Ax = Aoperator(Fx,Fy,uv,N,smoothness)
%     u = uv(1:N);
%     v = uv((N+1):end);
%     temp1 = 2*Fx.*Fy;
%     Bv = temp1.*v;
%     Bu = temp1.*u;
%     Au = smoothness*2*sobel(sobel(u,2),2,1) + smoothness*2*sobel(sobel(u,1),1,1) + Bu;
%     Av = smoothness*2*sobel(sobel(v,2),2,1) + smoothness*2*sobel(sobel(v,1),1,1) + Bv;
%     Ax = [Au + Bv; Bu + Av];
% end