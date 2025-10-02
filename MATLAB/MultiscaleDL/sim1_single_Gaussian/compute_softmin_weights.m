function [W_all, Rsum] = compute_softmin_weights(X, K, J, tau, padMode, normalize)
% X: [1 N KJ T] spatial
if nargin < 5, padMode = 'zero'; end
if nargin < 6, normalize = false; end
Xs = squeeze(X);         % [N x KJ x T]
[Ns, KJ, Tt] = size(Xs);
if KJ ~= K*J, error('KJ mismatch'); end
% reshape to [K J N T]
Xtemp = zeros(K, J, Ns, Tt);
for k = 1:K
    i1 = 1 + (k-1)*J; i2 = J + (k-1)*J;
    block = permute(Xs(:, i1:i2, :), [2 1 3]); % [J x N x T]
    Xtemp(k,:,:,:) = block;
end
Sshifts = [0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
L = size(Sshifts,1);
W_all = zeros(K, L, Ns, Tt-1);
Rsum = 0;
epsNorm = 1e-12;
for k = 1:K
    for t = 1:(Tt-1)
        Xt  = squeeze(Xtemp(k,:,:,t));    % [J x N]
        Xtp = squeeze(Xtemp(k,:,:,t+1));  % [J x N]
        rs = zeros(L,1);
        for sidx = 1:L
            s = Sshifts(sidx,:);
            Xs = shift2d(Xt, s(1), s(2), padMode);
            diff = (Xtp - Xs);
            rs(sidx) = sum(diff(:).^2);
            if normalize
                denom = sum(Xtp(:).^2) + epsNorm;
                rs(sidx) = rs(sidx) / denom;
            end
        end
        scaled = -rs ./ tau;
        m = max(scaled);
        ws = exp(scaled - m);
        ws = ws / sum(ws);
        % store scalar weights (broadcast over pixels)
        for sidx = 1:L
            W_all(k, sidx, :, t) = ws(sidx); % same for all pixels
        end
        lse = m + log(sum(exp(scaled - m)));
        r_soft = -tau * lse;
        Rsum = Rsum + r_soft;
    end
end
end
