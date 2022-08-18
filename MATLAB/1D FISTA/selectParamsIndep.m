function select_ind = selectParamsIndep(mse_indep,l1_norm)

[M,T] = size(mse_indep);
select_ind = zeros(T,1);
for time = 1:T
    crit = abs(l1_norm(:,time)*0.5).^2 + abs(mse_indep(:,time)).^2;
    select_ind(time) = find( (crit == min(crit)),1 );
end

end

