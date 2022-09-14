function awmv = computeAWMV_1D(x,var_theta)
%computeAWMV
% Inputs
% x- (z,x,az_var,rad_var) array of fitted coefficients
% var_theta- (az_var) array of azimuthal variances in dictionary

% Outputs
% awmv_az- amplitude weighted mean variance (azimuthal)

K = numel(var_theta);

% row vector
var_theta = var_theta(:)';

az_signal = squeeze(sum(x,1));
[n1,n2] = size(az_signal);
if n2 == K
    az_signal = az_signal';
end
total = sum(az_signal,1);




awmv = sqrt(var_theta)*az_signal./total;
