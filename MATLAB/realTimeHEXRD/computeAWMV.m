function [awmv_rad, awmv_az] = computeAWMV(x,P)
%computeAWMV
% Inputs
% x- (z,x,az_var,rad_var) array of fitted coefficients
% var_theta- (az_var) array of azimuthal variances in dictionary
% var_rad- (rad_var) array of radial variances in dictionary
% Outputs
% awmv_az- amplitude weighted mean variance (azimuthal)
% awmv_rad- amplitude weighted mean variance (rad)

% Sum over space, variance
x = reshape(x,[P.N1,P.N2,P.K1,P.K2]);
az_signal = squeeze(sum(x,[1,2,3]));
rad_signal = squeeze(sum(x,[1,2,4]));
total = sum(az_signal,'all');

% Compute awmv_rad
awmv_rad = P.sigma1(:)'*rad_signal./total;


% Compute awmv_az
awmv_az = P.sigma2(:)'*az_signal./total;



end