function Xcomp = filterbankVis(x,X,bank)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[hLo,hHi] = filters(bank);
L = size(X,2);
layers = log2(L);
Nx = numel(x);
Xcomp = zeros(Nx,L);

for j = 1:L
    y = X(:,j);
    binStr = dec2bin(j,layers);
    for i = 1:layers
        if str2num(binStr(i))
            y = applyFilter(y,hHi);
        else
            y = applyFilter(y,hLo);
        end
    end
    Ny = numel(y);
    diffSize = Ny-Nx;
    Xcomp(:,j) = y(1+floor(diffSize/2):end-ceil(diffSize/2));
end

end
function lo = applyFilter(in,hFilt)
    lo = conv(upsample(in,2),hFilt(:,2),'valid');
end
function lo = onlyFilter(in,hFilt)
    lo = conv(in,hFilt(:,2));
end


