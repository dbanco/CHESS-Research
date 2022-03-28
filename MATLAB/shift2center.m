function y = shift2center(x,p)
if nargin == 2
   if numel(size(x)) == 3
        [~,~,n3] = size(x);
        y = x;
        for i = 1:n3
            y(:,:,i) = shift2D(x(:,:,i),p(1),p(2));
        end
    elseif numel(size(x)) == 2
        y =shift2D(x,p(1),p(2));
    end 
else
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

end

