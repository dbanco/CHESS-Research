function D = createDiffMatrix(N,M,T)
%createDiffMatrix Summary of this function goes here
%   Detailed explanation goes here
D = zeros(T-1,T*N*M);
% row1 = [-ones(1,N*M), ones(1,N*M), zeros(1,(T-2)*N*M)];
for i = 1:(T-1)
    D(i,(1+(i-1)*N*M):(i)*N*M) = -1; 
    D(i,(1+(i)*N*M):(i+1)*N*M) =  1; 
end

