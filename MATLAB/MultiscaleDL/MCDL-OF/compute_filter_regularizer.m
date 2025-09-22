function [forwardOp] = compute_filter_regularizer(kernel,X,K,J,transpose)
% Assume
% X - [1,N,KJ,T]

[~,N,KJ,T] = size(X);

if nargin < 5
    transpose = false;
end

% Set up transpose of filter
if transpose
    padPattern = [1 1 1];
    for i = 1:numel(size(kernel))
        kernel = flip(kernel,i);
    end
    outT = T+1;
else 
    padPattern = [1 1];
    outT = T-1;
end

try
    forwardOp = gpuArray(complex(zeros(N,J,outT)));
catch
    forwardOp = zeros(N,J,outT);
end

for k = 1:K
    i1 = 1 + (k-1)*J;
    i2 = J + (k-1)*J;
    Xblock = reshape(X(1,:,i1:i2,:),[N,J,T]);
    Xpad = padarray(Xblock, padPattern, 0, 'both');
    forwardOp(:,i1:i2,:) = convn(Xpad,kernel,'valid');        
end

forwardOp = reshape(forwardOp,[1,N,KJ,outT]);

end