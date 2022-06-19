function [y,N,T] = two_shape_problem
%% Construct 1D test problem Box and Gaussian
T = 30;
N = 64;

boxPos = 15;
gausPos = 45;

widthT = [linspace(6,16,15),linspace(16,6,15)];
stdT = [linspace(3,10,15),linspace(10,3,15)];

y = zeros(N,T);
for t = 1:T
    y(:,t) = y(:,t) + makeBox(widthT(t),boxPos,N) +...
                      gaussian_basis_1D(N,gausPos,stdT(t));
end

end
function box = makeBox(width,pos,T)

box = zeros(T,1);
b1 = round(pos-width/2);
b2 = round(pos+width/2);
box(b1:b2) = 1;

end