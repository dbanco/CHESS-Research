function [awmv_az, awmv_rad] = computeAWMV(x,P)
%computeAWMV
% Inputs
% x- (z,x,az_var,rad_var) array of fitted coefficients
% var_theta- (az_var) array of azimuthal variances in dictionary
% var_rad- (rad_var) array of radial variances in dictionary
% Outputs
% awmv_az- amplitude weighted mean variance (azimuthal)
% awmv_rad- amplitude weighted mean variance (rad)

% Sum over space, variance
az_signal = squeeze(sum(sum(sum(x,1),2),4));
rad_signal = squeeze(sum(sum(sum(x,1),2),3));
total = sum(az_signal,'all');

% Compute awmv_rad
awmv_rad = P.sigma1'.*rad_signal./total;
awmv_rad = sum(awmv_rad(:));

% Compute awmv_az
awmv_az = P.sigma2'.*az_signal./total;
awmv_az = sum(awmv_az(:));


end