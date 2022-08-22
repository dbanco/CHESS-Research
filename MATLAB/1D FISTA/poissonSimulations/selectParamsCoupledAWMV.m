function select_ind = selectParamsCoupledAWMV(awmv,theta_stds)

[~,MM] = size(awmv);
awmv_err = zeros(MM,1);
for i = 1:MM
    awmv_err(i) = norm(awmv(:,i)-theta_stds(:));
end
select_ind = find( (awmv_err == min(awmv_err)),1 );

end

