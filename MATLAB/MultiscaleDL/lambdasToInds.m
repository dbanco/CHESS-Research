function selected_inds = lambdasToInds(lambdas,lamArrays)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = numel(lamArrays);
selected_inds = zeros(N,1);
for i = 1:N
    selected_inds(i) = find(lambdas(i)==lamArrays{i});
end


end