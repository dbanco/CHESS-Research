function DiffPhi = DiffPhiX_1D(X)
%DiffPhiX_1D Apply summing over space and temporal difference matrix 

[~,K,T] = size(X);
DiffPhi = zeros(K,T-1);
PhiX = squeeze(sum(X,1));
for t = 1:(T-1)
    DiffPhi(:,t) = PhiX(:,t+1)-PhiX(:,t);
end

