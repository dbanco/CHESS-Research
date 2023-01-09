function D = combineDictD(Dh,Dv,K)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

M1 = size(Dv,1);
M2 = size(Dh,2);
D = zeros(M1,M2,K);

for k = 1:K
     D(:,:,k) = Dv(:,:,k)*Dh(:,:,k);
end

end

