function evar = compute_exp_az_variance(x,var_theta)
%compute_exp_az_variance
% Inputs
% x- (z,x,az_var,rad_var) array of fitted coefficients
% var_theta- (az_var) array of azimuthal variances
% Outputs
% evar- expected variance scalar

% Sum over space, radial variance
az_signal = squeeze(sum(sum(sum(x,1),2),4));
total = sum(az_signal);

% Compute expected variance
evar = var_theta'.*az_signal./total;
evar = sum(evar(:));
end