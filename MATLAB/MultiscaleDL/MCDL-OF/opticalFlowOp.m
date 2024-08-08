function out = opticalFlowOp(x,u,v,K,trans)
%opticalFlowOp

if nargin < 5
    trans = 0;
end
if numel(size(x)) == 4
    reshapeIt = 1;
else
    reshapeIt = 0;
end

x = squeeze(x);
U = size(x,2)/K;

[Fx,Fy,Ft] = diffDict(x,K,U);

% Option to additionally apply operation transposed or just apply operation
% transposed
if trans == 1
    Ax = Fx.*u + Fy.*v + Ft;
    [uFx,vFy,Ft] = diffDictTrans(Ax,K,U,u,v);
    out = uFx + vFy + Ft;
elseif trans == 0
    out = Fx.*u + Fy.*v + Ft;
elseif trans == 2
    [uFx,vFy,Ft] = diffDictTrans(x,K,U,u,v);
    out = uFx + vFy + Ft;
end
if reshapeIt
    out = reshape(out,[1 size(out)]);
end

end