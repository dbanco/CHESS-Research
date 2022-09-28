function select_ind = selectParamsCoupled(mse,l1_norm,tv_penalty)

crit1 = 1.5*abs(mse-min(mse)).^2 /(max(mse)-min(mse))^2;
crit2 = abs(l1_norm-min(l1_norm)).^2 /2/(max(l1_norm)-min(l1_norm))^2;
crit3 = abs(tv_penalty-min(tv_penalty)).^2 /(max(tv_penalty)-min(tv_penalty))^2;
crit = crit1+crit2+crit3;
select_ind = find( (crit == min(crit)),1 );

end

