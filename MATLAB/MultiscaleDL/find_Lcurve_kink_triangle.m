function [idx_kink, lambda_kink, distances] = find_Lcurve_kink_triangle(residuals, solutions, lambda)
% FIND_LCURVE_KINK_TRIANGLE Finds the kink of an L-curve using the triangle method.
%
% INPUTS:
%   residuals  - vector of residual norms ||Ax - b|| for different lambda values
%   solutions  - vector of solution norms ||x|| for different lambda values
%   lambda     - vector of corresponding lambda values
%
% OUTPUTS:
%   idx_kink     - index of lambda corresponding to the kink (original input order)
%   lambda_kink  - lambda value at the kink
%   distances    - distances of each point from the straight line (sorted order)
%
% Example:
%   [idx, lam] = find_Lcurve_kink_triangle(residuals, solutions, lambda);

    if length(residuals) ~= length(solutions) || length(lambda) ~= length(solutions)
        error('residuals, solutions, and lambda must have the same length.');
    end
    
    % Sort by lambda and keep original indices
    [lambda_sorted, sortIdx] = sort(lambda(:));
    residuals_sorted = residuals(sortIdx);
    solutions_sorted = solutions(sortIdx);
    
    % Convert to log-log space
    r = log(residuals_sorted(:));
    s = log(solutions_sorted(:));
    
    % First and last points
    p1 = [r(1), s(1)];
    pN = [r(end), s(end)];
    
    % Line vector
    lineVec = pN - p1;
    
    % Compute perpendicular distance for each point
    distances = zeros(size(r));
    for i = 1:length(r)
        p = [r(i), s(i)];
        % Cross product magnitude formula for point-to-line distance
        distances(i) = det([lineVec; p - p1]) / norm(lineVec);
    end
    
    % Ignore first and last
    distances([1 end]) = 0;
    
    % Find max distance point (sorted index)
    [~, idx_sorted] = min(distances);
    
    % Map back to original index
    idx_kink = sortIdx(idx_sorted);
    lambda_kink = lambda(idx_kink);
end
