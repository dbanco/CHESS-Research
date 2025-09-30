function [gradR, Ddiag] = softmin_grad_majorizer(X, tau)
% Compute gradient and diagonal majorizer for softmin temporal regularizer
%
% Inputs:
%   X   : array of coefficients (K x S x Sigma x T)
%   tau : softmin temperature
%
% Outputs:
%   gradR : same size as X, gradient of R(X)
%   Ddiag : same size as X, diagonal majorizer entries (positive)

[K,S,SIG,T] = size(X);

gradR = zeros(K,S,SIG,T);
Ddiag = zeros(K,S,SIG,T);

% loop over k,s,sig,t (can vectorize later)
for k = 1:K
  for sig = 1:SIG
    for t = 1:(T-1)
      for s = 1:S
        % Collect neighbor diffs z_ij
        zvals = [];
        for di = -1:1
          for dj = -1:1
            si = s+di; sj = sig+dj;
            if si>=1 && si<=S && sj>=1 && sj<=SIG
              zvals(end+1) = X(k,s,sig,t) - X(k,si,sj,t+1);
            end
          end
        end
        
        % Soft weights
        expvals = exp(-(zvals.^2)/tau);
        w = expvals / sum(expvals);
        
        % Gradient contribution wrt X(k,s,sig,t)
        gradR(k,s,sig,t) = (2/tau) * sum(w .* zvals);
        
        % Diagonal majorizer: positive curvature upper bound
        Ddiag(k,s,sig,t) = (2/tau) * sum(w);
      end
    end
  end
end

end
