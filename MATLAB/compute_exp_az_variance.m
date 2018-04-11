function evar = compute_exp_az_variance(x,var_theta)
az_signal = squeeze(sum(sum(sum(x,1),2),4));
total = sum(az_signal);
evar = sum(var_theta.*az_signal/total);

end