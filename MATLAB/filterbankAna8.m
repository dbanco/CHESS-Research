function Xout = filterbankAna8(x,bank,layers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[hLo,hHi] = filters(bank);
M = numel(x);
L = 2^layers;
Xold = x;

for i = 1:layers
    L = 2^i;
    
    Xnew  = zeros(ceil((M+9)/2),L);
    for j = 1:L
        jj = floor((j+1)/2);
        if j<=L/2
            Xnew(:,j) = applyFilter(Xold(:,jj),hLo);
        else
            Xnew(:,j) = applyFilter(Xold(:,jj),hHi);
        end
    end
    Xold = Xnew;
    M = size(Xold,1);
end
Xout = Xnew;
end
function lo = applyFilter(in,hLo)
    lo = downsample(conv(in,hLo(:,1)),2);
end


