function b = gaussian_basis_wrap_1D(N_x,mean_x,std_x,scaling)
%gaussian_basis_wrap_1D Generates gaussian peak function vector
% Inputs:
% N_x - vector length
% mean_x - mean of gaussian basis function
% std_x - standard deviation of gaussian basis function
% scaling - '2-norm' unit 2-norm scaling
%           '1-norm' unit 1-norm scaling
%           'max'    unit max scaling factor
%           'rms'    unit root-mean-square
%           default is standard Gaussian scaling factor

% Outputs:
% b - (N x 1) vector

% Compute theta distances with wrapping at boundaries
idx = 1:N_x;
wrapN = @(x, N) (1 + mod(x-1, N));
opposite = (idx(wrapN(floor(mean_x-N_x/2),N_x)) +... 
            idx(wrapN(ceil(mean_x-N_x/2),N_x)))/2;
if opposite == mean_x
    opposite = 0.5;
end
dist1 = abs(mean_x - idx);
dist2 = N_x/2 - abs(opposite - idx);
dist = min(dist1,dist2);
dist_sq_theta = dist.^2;    % num_theta length vector

% Compute values
b = exp(-dist_sq_theta/(2*std_x^2) )';
if nargin > 3
    switch scaling
        case '2-norm'
            b = b/norm(b(:));
        case '1-norm'
            b = b/sum(abs(b(:)));
        case 'max'
            b = b/max(b(:));
        case 'rms'
            b = b/sqrt( sum(b(:).^2)/N_x );
        otherwise
            b = b/(std_x*sqrt(2*pi));
    end
end

end

