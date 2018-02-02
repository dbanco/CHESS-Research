function neighbors = neighbors_coef(x_hat_array,i,j,N,M)
%neighbors Returns sum of neighboring coefficients of (i,j)

neighbors = zeros(size(x_hat_array{1,1}));
count = 0;
for k = [i-1,i+1]
    if k>=1 && k<=N
        neighbors = neighbors + x_hat_array{k,j};
        count = count + 1;
    end
end
for k = [j-1,j+1]
    if k>=1 && k<=M
        neighbors = neighbors + x_hat_array{i,k};
        count = count + 1;
    end
neighbors = neighbors/count;

end

