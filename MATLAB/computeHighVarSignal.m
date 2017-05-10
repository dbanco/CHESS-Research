function [ high_az_var, high_rad_var ] = computeHighVarSignal( x_hat, az_idx, rad_idx )
%computeHighVarSignal Computes contribution of high variance bases to fit
% with respective to provided cutoffs
%
% Inputs:
% x_hat - (m x n x t x r) coefficients
% az_idx - index specifying lowest variance to include in high_az_var
% rad_idx - index specifying lowest variance to include in high_rad_var
%
% Outputs:
% high_az_var - proportion contributed by bases with high azimuthal variance
% high_rad_var - proportion contributed by bases with high radial variance

var_signal = squeeze(sum(sum(x_hat,1),2));
total_signal = sum(var_signal(:));
high_az_var = squeeze(sum( x_hat(az_idx:end, :) ))/total_signal;
high_rad_var = squeeze(sum( x_hat(:,rad_idx:end) ))/total_signal;

end
