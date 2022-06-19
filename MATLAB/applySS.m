function Y = applySS(X)
[T,R] = size(X);
Y = zeros(T,R);
Y(1,:) = X(2,:);
Y(T,:) = X(T-1,:);
for i = 2:T-1
    Y(i,:) = 0.5*(X(i-1,:) + X(i+1,:));
end
end
