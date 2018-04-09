function evar = compute_exp_az_variance(var_signal,var_theta)
az_var_signal = squeeze(sum(var_signal,2));
[~, s2, s3, s4] = size(az_var_signal);
evar = zeros(s2,s3,s4);
for j = 1:s2
    for k = 1:s3
        for l = 1:s4
            total = sum(az_var_signal(:,j,k,l));
            for i = 1:numel(var_theta)
                evar(j,k,l) = evar(j,k,l) +...
                      var_theta(i)*az_var_signal(i,j,k,l)/total;
            end
        end        
    end
end