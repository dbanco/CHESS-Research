function select_ind = selectParamsIndep(mse_indep,l1_norm)

[T,M] = size(mse_indep);
select_ind = zeros(T,1);
for time = 1:T
    crit = 0.5*abs(l1_norm(time,:)).^2 + abs(mse_indep(time,:)).^2;
    select_ind(time) = find( (crit == min(crit)),1 );
end

end

