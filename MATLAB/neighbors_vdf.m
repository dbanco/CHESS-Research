function [neighbors, vdfs] = neighbors_vdf(x_hat_array,i,j)
%neighbors Returns term used in gradient of vdf objective and vdfs 
%          of coefficients neighboring point (i,j)
[N,M] = size(x_hat_array);
[n,m,T,R] = size(x_hat_array{1,1});
neighbors = zeros(n,m,T,R);
vdfs = {}; idx = 1;
count = 0;
% Vertical neighbors
for k = [i-1,i+1]
    if k>=1 && k<=N
        x_hat_var = x_hat_array{k,j};
        vdf = squeeze(sum(sum(x_hat_var)))/sum(x_hat_var(:));
        vdfs{idx} = vdf;
        idx = idx + 1;
        for t = 1:T
            for r = 1:R
                neighbors(:,:,t,r) = neighbors(:,:,t,r) + ones(n,m)*vdf(t,r);
            end
        end
    end
end
% Horizontal neighbors
for k = [j-1,j+1]
    if k>=1 && k<=M
        x_hat_var = x_hat_array{i,k};
        vdf = squeeze(sum(sum(x_hat_var)))/sum(x_hat_var(:));
        vdfs{idx} = vdf;
        idx = idx + 1;
        for t = 1:T
            for r = 1:R
                neighbors(:,:,t,r) = neighbors(:,:,t,r) + ones(n,m)*vdf(t,r);
            end
        end
    end

end

