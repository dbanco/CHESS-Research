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

wrapN = @(x, N) (1 + mod(x-1, N));

% Compute theta distances with wrapping at boundaries
idx = 1:N_x;

mean_x = wrapN(mean_x,N_x);

upDiff = ceil(mean_x) - mean_x;
downDiff = mean_x - floor(mean_x);

if (mod(mean_x,1)==0) & (upDiff == downDiff)
    upDist = idx - 1 + upDiff;
    downDist = flip( idx - 1 +downDiff);

    dist1 = circshift(upDist,floor(mean_x)-1);
    dist2 = circshift(downDist,floor(mean_x));   
else
    upDist = idx - 1 + upDiff;
    downDist = flip( idx - 1 +downDiff);

    dist1 = circshift(upDist,floor(mean_x));
    dist2 = circshift(downDist,floor(mean_x));
end

% Compute values
b1 = exp(-dist1.^2/(2*std_x^2) )';
b2 = exp(-dist2.^2/(2*std_x^2) )';
b = b1 + b2;
if floor(mean_x) == mean_x
    b(mean_x) = b1(mean_x);
end
for i = 1:floor(6*std_x/N_x)
    b = b + exp(-(dist1+i*N_x).^2/(2*std_x^2) )' +...
            exp(-(dist2+i*N_x).^2/(2*std_x^2) )';
end
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

