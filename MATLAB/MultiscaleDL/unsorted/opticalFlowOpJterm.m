function out = opticalFlowOpJterm(x,u,v,K,trans,Jterm)
%opticalFlowOp

if nargin < 4
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

% Jx = x;
% if nargin < 6
%     Jx(Ft==0) = 0;
% else
%     Jx(Jterm) = 0;
% end

% Option to additionally apply operation transposed
if trans
    Ax = Fx.*u + Fy.*v + Ft;
%     Ax = Fx.*u + Fy.*v + Ft + Jx;
    [Fx,Fy] = sobelDict(Ax,K,U,trans);
    Ft = timeDiffDict(Ax,K,U,trans);
    
%     Jx = Ax;
%     if nargin < 6 
%         Jx(Ft==0) = 0;
%     else
%         Jx(Jterm) = 0;
%     end
    out = Fx.*u + Fy.*v + Ft;
%     out = Fx.*u + Fy.*v + Ft + Jx;
else
    out = Fx.*u + Fy.*v + Ft;
%     out = Fx.*u + Fy.*v + Ft + Jx;
end
if reshapeIt
    out = reshape(out,[1 size(out)]);
end

end