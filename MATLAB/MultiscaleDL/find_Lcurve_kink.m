function [idx_kink, lambda_kink] = find_Lcurve_kink(residuals, solutions, lambda, smooth_win)
% FIND_LCURVE_KINK Finds the L-curve corner using max curvature after sorting and smoothing.
%
% INPUTS:
%   residuals   - vector of residual norms ||Ax - b|| for different lambda values
%   solutions   - vector of solution norms ||x|| for different lambda values
%   lambda      - vector of corresponding lambda values
%   smooth_win  - optional smoothing window size (default: 5)
%
% OUTPUT:
%   idx_kink    - index of lambda corresponding to maximum curvature (after sorting)
%   lambda_kink - lambda value at the kink
%
% Example:
%   [idx, lam] = find_Lcurve_kink(residuals, solutions, lambda, 7);

    if nargin < 4
        smooth_win = 5; % default smoothing window
    end

    if length(residuals) ~= length(solutions) || length(lambda) ~= length(solutions)
        error('residuals, solutions, and lambda must have the same length.');
    end
    
    % Sort all inputs by lambda (in ascending order)
    [lambda, sortIdx] = sort(lambda(:));
    residuals = residuals(sortIdx);
    solutions = solutions(sortIdx);
    
    % Convert to log scale
    r = log(residuals(:));
    s = log(solutions(:));
    
    % Smooth the data
    r = smooth(r, smooth_win);
    s = smooth(s, smooth_win);
    
    % Compute first and second derivatives
    dr = gradient(r);
    ds = gradient(s);
    d2r = gradient(dr);
    d2s = gradient(ds);
    
    % Compute curvature: |r'' * s' - r' * s''| / (r'^2 + s'^2)^(3/2)
    curvature = abs(d2r .* ds - dr .* d2s) ./ ((dr.^2 + ds.^2).^(3/2));
    
    % Find maximum curvature
    [~, idx_kink] = max(curvature);
    lambda_kink = lambda(idx_kink);
    idx_kink = sortIdx(idx_kink);
end
