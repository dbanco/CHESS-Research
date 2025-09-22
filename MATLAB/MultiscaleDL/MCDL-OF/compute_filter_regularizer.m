function [forwardOp] = compute_filter_regularizer(kernel,X,K,J,transpose)
% Assume
% X - [1,N,KJ,T]

if nargin < 5
    transpose = false;
end

if transpose
    % Flip each dimenension for transposed convolution
    for i = 1:numel(size(filter))
        kernel = flip(kernel,i);
    end
end

forwardOp = zeros(size(X));
for k = 1:K
    i1 = 1 + (k-1)*J;
    i2 = J + (k-1)*J;
    forwardOp(1,:,i1:i2,:) = convn(squeeze(X(1,:,i1:i2,:)),kernel,'same');        
end

end