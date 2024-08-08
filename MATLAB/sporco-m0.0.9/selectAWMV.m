function awmv_select = selectAWMV(awmv,select_ind);

[T,M] = size(awmv);
awmv_select = zeros(T,1);

for i = 1:T
    awmv_select(i) = awmv(i,select_ind(i));
end

end