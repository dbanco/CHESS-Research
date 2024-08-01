function [AG,NormsHV] = combineDict(N1,N2,AGh,AGv,K,U,NormsHV)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

M1 = size(AGv,1);
M2 = size(AGh,2);
AG = zeros(N1,N2,K,U);

AGv = reshape(AGv,[M1,1,K,U]);
AGh = reshape(AGh,[1,M2,K,U]);

for k = 1:K
    for u = 1:U
            AG(1:M1,1:M2,k,u) = AGv(:,:,k,u)*AGh(:,:,k,u);
    end
end
AG = reshape(AG,[N1,N2,K*U]);
if nargin < 7
    NormsHV = zeros(K*U,1);
    for i = 1:K*U
        NormsHV(i) = norm(AG(:,:,i),'fro');
    end
end
for i = 1:K*U
    AG(:,:,i) = AG(:,:,i)./NormsHV(i);
end

end

