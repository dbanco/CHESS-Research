function lam_inds = convertToInds(lam_vals,lam1,lam2,lam3)

N = size(lam_vals,1);
lam_inds = zeros(N,3);
for i = 1:N
    lam_inds(i,1) = find(lam_vals(i,1) == lam1);
    lam_inds(i,2) = find(lam_vals(i,2) == lam2);
    lam_inds(i,3) = find(lam_vals(i,3) == lam3);
end

end