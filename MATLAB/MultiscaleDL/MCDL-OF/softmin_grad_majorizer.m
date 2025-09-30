function [gradR_out, Ddiag_out] = softmin_grad_majorizer(X, K, J, tau)
% Compute gradient and diagonal majorizer for softmin temporal regularizer
% Numerically stable version with vectorized exponentials
%
% Inputs:
%   X   : array of size [1, N+M-1, K*J, T]
%   K   : number of atoms
%   J   : number of channels per atom
%   tau : softmin temperature (positive scalar)
%
% Outputs:
%   gradR : same size as reshaped X [K x J x N x T], gradient of R(X)
%   Ddiag : same size as X, diagonal majorizer entries (positive)

[~, N, KJ, T] = size(X);

% Reshape to [K x J x N x T]
Xtemp = zeros(K, J, N, T);
Xs = squeeze(X);  
for k = 1:K
    i1 = 1 + (k-1)*J;
    i2 = J + (k-1)*J;
    % Xs(:, i1:i2, :) is [N x J x T] -> permute to [J x N x T]
    block = permute(Xs(:, i1:i2, :), [2 1 3]); 
    Xtemp(k,:,:,:) = block;
end
X = Xtemp;

gradR = zeros(size(X));
Ddiag = zeros(size(X));

eps_small = 1e-12;  % prevent division by zero

% Loop over temporal slices
for t = 1:(T-1)
    % Current and next slice
    X_curr = X(:,:,:,t);       % [K x J x N]
    X_next = X(:,:,:,t+1);     % [K x J x N]
    
    % Pad next slice along spatial dims (J and N)
    X_next_pad = padarray(X_next, [0 1 1], 'replicate'); % pad J,N
    
    % Compute 3x3 neighborhood differences
    zvals = zeros(K, J, N, 9);  % 9 neighbors
    n = 1;
    for di = -1:1
        for dj = -1:1
            zvals(:,:,:,n) = X_curr - X_next_pad(:, (1:J)+di+1, (1:N)+dj+1);
            n = n + 1;
        end
    end
    
    % Numerically stable softmin weights
    z2_tau = -(zvals.^2) / tau;
    z2_tau = z2_tau - max(z2_tau, [], 4);   % subtract max along neighbor dim
    expvals = exp(z2_tau);
    w = expvals ./ (sum(expvals, 4) + eps_small);
    
    % Gradient
    gradR(:,:,:,t) = (2/tau) * sum(w .* zvals, 4);
    
    % Diagonal majorizer
    Ddiag(:,:,:,t) = (2/tau) * sum(w, 4);
end

% Last temporal slice: gradient = 0, majorizer = previous slice
gradR(:,:,:,T) = 0;
Ddiag(:,:,:,T) = Ddiag(:,:,:,T-1);

gradR_out = zeros(1,N,KJ,T);
Ddiag_out = zeros(1,N,KJ,T);
for k = 1:K
    i1 = 1 + (k-1)*J;
    i2 = J + (k-1)*J;
    block = permute(gradR(k,:,:,:), [1 3 2 4]); 
    gradR_out(1,:,i1:i2,:) = block;
    block = permute(Ddiag(k,:,:,:), [1 3 2 4]); 
    Ddiag_out(1,:,i1:i2,:) = block;
end

end
