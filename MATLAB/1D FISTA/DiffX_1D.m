function Diff = DiffX_1D(X)
%DiffX_1D Apply temporal difference matrix

[N,K,T] = size(X);
Diff = zeros(N,K,T-1);
for t = 1:(T-1)
    Diff(:,:,t) = X(:,:,t+1)-X(:,:,t);
end

