function b = voigt_basis_wrap_1D(N,mu,gamma,eta,scaling)
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
dist_sq = dist.^2;    % num_theta length vector

sigma = gamma/(2*sqrt(2*log(2)));
ag = 1./(sigma*sqrt(2*pi));
bg = 4*log(2)/gamma^2;

% Compute values
b_lorentzian = (gamma/(2*pi))./(dist_sq + (gamma/2)^2 )';
b_gaussian = ag*exp(-bg*dist_sq/(sigma^2) )';
b = eta*b_gaussian + (1-eta)*b_lorentzian;

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
    end
end

end

