function [awmv_az, awmv_rad] = computeAWMV(x,var_theta,var_rad)
%computeAWMV
% Inputs
% x- (z,x,az_var,rad_var) array of fitted coefficients
% var_theta- (az_var) array of azimuthal variances in dictionary
% var_rad- (rad_var) array of radial variances in dictionary
% Outputs
% awmv_az- amplitude weighted mean variance (azimuthal)
% awmv_rad- amplitude weighted mean variance (rad)

% Sum over space, radial variance
az_signal = squeeze(sum(sum(sum(x,1),2),4));
rad_signal = squeeze(sum(sum(sum(x,1),2),3));
total = sum(az_signal);

% Compute aemv_az
awmv_az = var_theta'.*az_signal./total;
awmv_az = sum(awmv_az(:));

% Compute aemv_rad
awmv_rad = var_rad'.*rad_signal./total;
awmv_rad = sum(awmv_rad(:));
end