function Y = applyL(X)
[T,R] = size(X);
Y = zeros(T-1,R);
for t = 1:(T-1)
    Y(t,:) = X(t,:) - X(t+1,:);
end
end
