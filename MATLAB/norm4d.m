function y = norm4d( x, p )
%norm4d Summary of this function goes here
%   Detailed explanation goes here
[a,b,c,d] = size(x);
y = norm(reshape(x,a*b,c*d),p);
end

