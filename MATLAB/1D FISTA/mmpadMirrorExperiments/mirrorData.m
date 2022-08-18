function b_m = mirrorData(b)
nn = numel(b);
pad1 = floor(nn/2);
pad2 = ceil(nn/2);
N = nn + pad1 + pad2;
b_m = zeros(N,1);
b_m((pad1+1):(pad1+nn)) = b;
b_m((1+N-pad2):N) = flip(b((nn-pad2+1):nn));
b_m(1:pad1) = flip(b(1:pad1));
end

