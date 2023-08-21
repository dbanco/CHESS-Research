function x = unpad(xpad,M,shape)
if numel(size(xpad)) == 2
    switch shape
        case 'both'
            x = xpad(1+M:end-M,:);
        case 'post'
            x = xpad(1:end-M,:);
        case 'pre'
            x = xpad(1+M:end,:);
    end
elseif numel(size(xpad)) == 4
    switch shape
        case 'both'
            x = xpad(:,1+M:end-M,:,:);
        case 'post'
            x = xpad(:,1:end-M,:,:);
        case 'pre'
            x = xpad(:,1+M:end,:,:);
    end
end
end