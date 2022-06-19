function D = initGaussianDictionary(dictFilterSizes)
K = size(dictFilterSizes,2);
P.N = max(dictFilterSizes(:));
P.stds = logspace(log10(0.5),log10(P.N/6),K);
P.means = dictFilterSizes(2,:)/2;
P.basis = 'norm2';
D = zeros(1,P.N,K);
for i = 1:K
    D(1,:,:) = dictionary(P);
end
end