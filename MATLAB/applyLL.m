function Y = applyLL(X)
[T,R] = size(X);
Y = zeros(T,R);
Y(1,:) = -2*X(1,:) + 2*X(2,:);
Y(T,:) = 2*X(T-1,:) - 2*X(T,:);
for t = 2:(T-1)
    Y(t,:) = -2*X(t,:) + X(t+1,:) + X(t-1,:);
end
end
