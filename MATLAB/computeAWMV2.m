function [awmv_az, awmv_rad] = computeAWMV2(x,var_theta,var_rad,norms)
%computeAWMV
% Inputs
% x- (z,x,az_var,rad_var) array of fitted coefficients
% var_theta- (az_var) array of azimuthal variances
% Outputs
% evar- expected variance scalar

% Sum over space, radial variance
az_signal = squeeze(sum(sum(sum(x,1),2),4));
rad_signal = squeeze(sum(sum(sum(x,1),2),3));
total = sum(az_signal);

% Compute aemv_az
awmv_az = norms'*(var_theta.*az_signal./total);
awmv_az = sum(awmv_az(:));

% Compute aemv_rad
awmv_rad = norms*(rad_signal.*var_rad./total);
awmv_rad = sum(awmv_rad(:));
end