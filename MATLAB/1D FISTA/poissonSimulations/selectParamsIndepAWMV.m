function select_ind = selectParamsIndepAWMV(awmv,theta_stds)

[T,M] = size(awmv);
select_ind = zeros(T,1);
awmv_err = zeros(T,M);
for i = 1:M
    awmv_err(:,i) = abs(awmv(:,i)-theta_stds(:));
end
for t = 1:T
    amwv_err_t = awmv_err(t,:);
    select_ind(t) = find( (amwv_err_t == min(amwv_err_t)),1 );
end

