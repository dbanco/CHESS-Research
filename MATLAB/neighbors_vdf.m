function neighbors = neighbors_vdf(x_hat_array,i,j)
%neighbors Returns sum of neighboring coefficients of (i,j)
[N,M] = size(x_hat_array);
[n,m,T,R] = size(x_hat_array{1,1});
neighbors = zeros(n,m,T,R);
count = 0;
for k = [i-1,i+1]
    if k>=1 && k<=N
        x_hat_var = x_hat_array{k,j};
        vdf = squeeze(sum(sum(x_hat_var)))/sum(x_hat_var(:));
        for t = 1:T
            for r = 1:R
                neighbors(:,:,t,r) = neighbors(:,:,t,r) + ones(n,m)*vdf(t,r);
            end
        end
    end
end
for k = [j-1,j+1]
    if k>=1 && k<=M
        x_hat_var = x_hat_array{i,k};
        vdf = squeeze(sum(sum(x_hat_var)))/sum(x_hat_var(:));
        for t = 1:T
            for r = 1:R
                neighbors(:,:,t,r) = neighbors(:,:,t,r) + ones(n,m)*vdf(t,r);
            end
        end
    end

end

