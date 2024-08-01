function out = opticalFlowOpHS(x,u,v,K,trans)
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

[Fx,Fy] = sobelDict(x,K,U);
Ft = timeDiffDict(x,K,U);

% Option to additionally apply operation transposed
if trans
    Ax = Fx.*u + Fy.*v + Ft;
    [uFx,vFy] = sobelDictTrans(Ax,K,U,u,v);
    Ft = timeDiffDict(Ax,K,U,trans);
    out = uFx + vFy + Ft;
else
    out = Fx.*u + Fy.*v + Ft;
end
if reshapeIt
    out = reshape(out,[1 size(out)]);
end

end