function awmv = computeAWMV_1D(x,var_theta)
%computeAWMV
% Inputs
% x- (z,x,az_var,rad_var) array of fitted coefficients
% var_theta- (az_var) array of azimuthal variances in dictionary

% Outputs
% awmv_az- amplitude weighted mean variance (azimuthal)


% row vector
var_theta = var_theta(:)';

az_signal = squeeze(sum(x,1));
total = sum(az_signal,1);
awmv = sqrt(var_theta)*az_signal./total;
