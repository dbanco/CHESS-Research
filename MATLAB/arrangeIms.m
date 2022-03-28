function dispX = arrangeIms(X)

if numel(size(X)) == 3
    [n1,n2,n3] = size(X);
    dispX = zeros(n1*n3,n2);
    for t = 1:n3
        dispX(1+n1*(t-1):t*n1,:) = X(:,:,t)/max(max(X(:,:,t)));
    end
elseif numel(size(X)) ==2
    [n1,n2] = size(X);
    dispX = zeros(n1,n2);
    for t = 1:n2
        dispX(:,t) = X(:,t)/max(max(X(:,t)));
    end    
end


end

