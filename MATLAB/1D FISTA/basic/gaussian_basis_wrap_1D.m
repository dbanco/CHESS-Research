function b = gaussian_basis_wrap_1D(N,mu,sigma,scaling)
%gaussian_basis_wrap_1D Generates gaussian peak function vector
% Inputs:
% N - vector length
% mu - mean of gaussian basis function
% sigma - standard deviation of gaussian basis function
% scaling - '2-norm' unit 2-norm scaling
%           '1-norm' unit 1-norm scaling
%           'max'    unit max scaling factor
%           'rms'    unit root-mean-square
%           default is standard Gaussian scaling factor

% Outputs:
% b - (N x 1) vector

% Compute theta distances with wrapping at boundaries
idx = 1:N;
wrapN = @(x, N) (1 + mod(x-1, N));
opposite = wrapN((mu-N/2),N);
            
if opposite == mu
    opposite = 0.5;
end
dist1 = abs(mu - idx);
dist2 = abs(N/2 - abs(opposite - idx));
dist = min(dist1,dist2);
dist_sq_theta = dist.^2;    % num_theta length vector

% Compute values
b = exp(-dist_sq_theta/(sigma^2) )';
if nargin > 3
    switch scaling
        case '2-norm'
            b = b/norm(b(:));
        case '1-norm'
            b = b/sum(abs(b(:)));
        case 'max'
            b = b/max(b(:));
        case 'rms'
            b = b/sqrt( sum(b(:).^2)/N );
        otherwise
            b = b/(sigma*sqrt(2*pi));
    end
end

end

