function [u,v,Fy,Fx,Ft] = computeHornSchunkLS(data,u0,v0,smoothness,maxIters,tol)
%computeHornSchunk Summary of this function goes here
%   Detailed explanation goes here

if nargin <4 || isempty(smoothness)
    smoothness = 1;
end
if nargin <5 || isempty(maxIters)
    maxIters = 100;
end
if nargin < 6 || isempty(tol)
  tol = 1e-4;
end

% Compute partial derivatives
Fy = sobel(data,1);
Fx = sobel(data,2);
Ft = -timeDiff(data);


% Set up linear system
b = [-2*Fx(:).*Ft(:); -2*Fy(:).*Ft(:)];
N = numel(data);
Aop = @(uv) Aoperator(Fx(:),Fy(:),uv,N,smoothness);

% Solve with conjugate gradient
uv0 = [u0(:);v0(:)];
[uv,~] = pcg(@(u) Aop(u), b(:),tol,maxIters,[],[],uv0);
u = reshape(uv(1:N),    size(data));
v = reshape(uv(N+1:end),size(data));

end

function Ax = Aoperator(Fx,Fy,uv,N,smoothness)
    u = uv(1:N);
    v = uv((N+1):end);
    temp1 = 2*Fx.*Fy;
    Bv = temp1.*v;
    Bu = temp1.*u;
    Au = smoothness*2*(diffx(diffx(u),1) + diffy(diffy(u),1)) + Bu;
    Av = smoothness*2*(diffx(diffx(v),1) + diffy(diffy(v),1)) + Bv;
    Ax = [Au + Bv; Bu + Av];
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