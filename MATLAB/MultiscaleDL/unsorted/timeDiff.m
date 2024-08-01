function dt = timeDiff(x,trans)
%timeDiff Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    trans = 0;
end
% append = zeros(size(x(:,:,1)));
append = x(:,:,1);
% append = x(:,:,end);
if trans == 0
    dt = (x(:,:,2:end)-x(:,:,1:end-1));
    dt = cat(3, append, dt);
elseif trans == 1
    dt = x;
    dt(:,:,1:end-1) = (dt(:,:,1:end-1) - x(:,:,2:end));
end