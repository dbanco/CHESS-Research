function [u,v,Fy,Fx,Ft] = computeLucasKanade(data,winSize)
%computeHornSchunk Summary of this function goes here
%   Detailed explanation goes here

if nargin <2
    winSize = 1;
end
[Ny,Nx,Nt] = size(data);

Fy = sobel(data,1);
Fx = sobel(data,2);
Ft = timeDiff(data);

u = zeros(size(data));
v = zeros(size(data));

for t = 1:Nt
    for i = 1:Ny
        for j = 1:Nx
            yWin = i-winSize:i+winSize;
            xWin = j-winSize:j+winSize;
            if xWin(1) < 1
                xWin = xWin + 1;
            end
            if yWin(1) < 1
                yWin = yWin + 1;
            end
            if yWin(end) > Ny
                yWin = yWin - 1;
            end
            if xWin(end) > Nx
                xWin = xWin - 1;
            end
            A = [vec(Fx(yWin,xWin,t)),vec(Fy(yWin,xWin,t))];
            b = - vec(Ft(yWin,xWin,t));
            if sum(abs(b)) > 100e-16
                uv = linsolve(A,b);
            else
                uv = [0,0];
            end
            u(i,j,t) = uv(1);
            v(i,j,t) = uv(2);
        end
    end
end

% constraint = sum(u.*Fx + v.*Fy + Ft,3);
% imagesc(constraint)
% norm(constraint(:))
% u = u;
% v = v;

end