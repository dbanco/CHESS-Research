function term = HornSchunckTerm(u,v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
term = sum(diffx(u(:)).^2,'all') + sum(diffy(u(:)).^2,'all') +...
       sum(diffx(v(:)).^2,'all') + sum(diffy(v(:)).^2,'all');
end