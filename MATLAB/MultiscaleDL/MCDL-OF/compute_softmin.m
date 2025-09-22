function [R, dR] = compute_softmin(X, K, J, tau, padMode, normalize)
% compute_softmin
% Soft-min temporal shift penalty
%
% Inputs
%   X       : [1 x N x (K*J) x T] coefficient tensor (double)
%   K       : number of atoms
%   J       : number of scales (so K*J is the 3rd dim block size)
%   tau     : soft-min temperature (suggest 1e-3 .. 1e-1). Default 1e-3.
%   padMode : 'zero' or 'reflect' (default 'zero')
%   normalize : boolean, if true normalizes each r_s by (||Xtp||^2 + eps). Default false.
%
% Outputs
%   R       : scalar penalty (sum over k and t)
%   dR_flat : gradient in same shape as X_flat (empty if not requested)
%
% Example:
%   [R, dR] = compute_time_reg_softmin_flat(X_flat, K, J, 1e-3, 'zero', true);

if nargin < 4 || isempty(tau), tau = 1e-3; end
if nargin < 5 || isempty(padMode), padMode = 'zero'; end
if nargin < 6 || isempty(normalize), normalize = false; end

% squeeze singleton leading dim
Xs = squeeze(X);         % [N x (K*J) x T]
[N, KJ, T] = size(Xs);
if KJ ~= K * J
    error('Third dim mismatch: expected K*J columns.');
end

% reshape to [K x J x N x T]
Xtemp = zeros(K, J, N, T);
for k = 1:K
    i1 = 1 + (k-1)*J;
    i2 = J + (k-1)*J;
    % Xs(:, i1:i2, :) is [N x J x T] -> permute to [J x N x T]
    block = permute(Xs(:, i1:i2, :), [2 1 3]); 
    Xtemp(k,:,:,:) = block;
end

% compute soft-min on [K x J x N x T] core
[R_core, dR_core] = compute_softmin_core(Xtemp, tau, padMode, normalize);

% map gradient back to flattened shape if requested
if nargout > 1 && ~isempty(dR_core)
    dR = zeros(1, N, KJ, T);
    for k = 1:K
        i1 = 1 + (k-1)*J;
        i2 = J + (k-1)*J;
        block = squeeze(dR_core(k,:,:,:)); % [J x N x T]
        block = permute(block, [2 1 3]);    % [N x J x T]
        dR(1, :, i1:i2, :) = block;
    end
else
    dR = [];
end

R = R_core;

end


%% ---------- CORE (works on [K x J x N x T]) ----------
function [R, dR] = compute_softmin_core(X, tau, padMode, normalize)
% X : [K x J x N x T]

if nargin < 2 || isempty(tau), tau = 1e-3; end
if nargin < 3 || isempty(padMode), padMode = 'zero'; end
if nargin < 4 || isempty(normalize), normalize = false; end

[K, J, N, T] = size(X);

% 8-neighborhood + identity
S = [ 0  0;
      1  0;
     -1  0;
      0  1;
      0 -1;
      1  1;
      1 -1;
     -1  1;
     -1 -1];

L = size(S,1);
R = 0;
computeGrad = nargout > 1;
if computeGrad, dR = zeros(size(X)); else dR = []; end

epsNorm = 1e-12;

% precompute j/n grids for any shift op if using vectorized shift helper (we use shift2d)
for k = 1:K
    for t = 1:(T-1)
        Xt  = squeeze(X(k,:,:,t));    % [J x N]
        Xtp = squeeze(X(k,:,:,t+1));  % [J x N]

        rs = zeros(L,1);
        Xshifts = cell(L,1);
        for sidx = 1:L
            s = S(sidx,:);
            Xs = shift2d(Xt, s(1), s(2), padMode); % shift Xt
            Xshifts{sidx} = Xs;
            diff = (Xtp - Xs);
            rs(sidx) = sum(diff(:).^2);
            if normalize
                denom = sum(Xtp(:).^2) + epsNorm;
                rs(sidx) = rs(sidx) / denom;
            end
        end

        % stable soft-min r_soft = -tau * logsumexp(-rs/tau)
        scaled = -rs ./ tau;
        m = max(scaled);
        lse = m + log(sum(exp(scaled - m)));
        r_soft = -tau * lse;
        R = R + r_soft;

        if computeGrad
            % weights (softmax over -rs/tau)
            ws = exp(scaled - m);
            ws = ws / sum(ws);    % Lx1

            grad_Xtp = zeros(J,N);
            grad_Xt  = zeros(J,N);

            for sidx = 1:L
                w = ws(sidx);
                Xs = Xshifts{sidx};
                diff = (Xtp - Xs);            % derivative wrt Xtp: 2*diff
                if normalize
                    denom = sum(Xtp(:).^2) + epsNorm;
                    % careful: rs used normalized; derivative includes denom terms
                    % r_s = sum(diff.^2)/denom
                    % dr/dXtp = (2*diff*denom - sum(diff.^2)*2*Xtp) / denom^2
                    ssum = sum(diff(:).^2);
                    grad_r_wrt_Xtp = (2*diff*denom - ssum*2*Xtp) / (denom^2);
                    grad_r_wrt_Xs = - (2*diff) / denom; % derivative wrt Xs (before shift-back)
                else
                    grad_r_wrt_Xtp = 2 * diff;
                    grad_r_wrt_Xs  = -2 * diff;
                end

                grad_Xtp = grad_Xtp + w .* grad_r_wrt_Xtp;

                % gradient wrt Xt: grad_r_wrt_Xs shifted back by -s
                s = S(sidx,:);
                grad_unshift = shift2d(grad_r_wrt_Xs, -s(1), -s(2), padMode);
                grad_Xt = grad_Xt + w .* grad_unshift;
            end

            % accumulate into full gradient tensor
            dR(k,:,:,t)   = squeeze(dR(k,:,:,t))   + grad_Xt;
            dR(k,:,:,t+1) = squeeze(dR(k,:,:,t+1)) + grad_Xtp;
        end
    end
end

if ~computeGrad
    dR = [];
end
end


%% ---------- helper: shift2d ----------
function Y = shift2d(X, shift_j, shift_n, padMode)
% shift2d shifts matrix X by (shift_j, shift_n) with padding
[J, N] = size(X);
Y = zeros(J, N);

j_src = (1:J) - shift_j;
n_src = (1:N) - shift_n;

valid_j = (j_src >= 1) & (j_src <= J);
valid_n = (n_src >= 1) & (n_src <= N);

if strcmp(padMode,'reflect')
    j_src = reflect_index(j_src, J);
    n_src = reflect_index(n_src, N);
    % vectorized assignment
    Y(:,:) = X(j_src, n_src);
else
    % zero padding
    for jj = 1:J
        if ~valid_j(jj), continue; end
        rs = j_src(jj);
        for nn = 1:N
            if ~valid_n(nn), continue; end
            cs = n_src(nn);
            Y(jj,nn) = X(rs, cs);
        end
    end
end
end

function idxr = reflect_index(idx, L)
% symmetric reflection into 1..L for vectorized mapping
idxr = idx;
idxr(idx < 1) = 2 - idx(idx < 1);
idxr(idx > L) = 2*L - idx(idx > L);
idxr = max(1, min(L, idxr));
end
