function [nablaU2,nablaV2] = laplaceOp(u,v,K,U)

avgKernel = [1 2 1; 2 0 2; 1 2 1]/12;

% If there are multiple atoms
if nargin > 2
    nablaU2 = zeros(size(u));
    nablaV2 = zeros(size(v));
    
    for k = 1:K
        i1 = 1+(k-1)*U;
        i2 = k*U;
        uk = u(:,i1:i2,:);
        vk = v(:,i1:i2,:);
        uAvg= convn(uk,avgKernel,'same');
        vAvg = convn(vk,avgKernel,'same');
        
        nablaU2(:,i1:i2,:) = uAvg - u(:,i1:i2,:);
        nablaV2(:,i1:i2,:) = vAvg - v(:,i1:i2,:);
    end
else
    uAvg = convn(u,avgKernel,'same');
    vAvg = convn(v,avgKernel,'same');
    
    nablaU2 = uAvg - u;
    nablaV2 = vAvg - v;
end


end