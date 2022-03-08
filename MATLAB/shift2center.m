function y = shift2center(x)
if numel(size(x)) == 3
    [n1,n2,n3] = size(x);
    y = x;
    for i = 1:n3
        y(:,:,i) = shift2D(x(:,:,i),floor(n1/2),floor(n2/2));
    end
elseif numel(size(x)) == 2
    [n1,n2] = size(x);
    y =shift2D(x,floor(n1/2),floor(n2/2));
end


end

