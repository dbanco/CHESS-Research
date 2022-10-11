function f = checkDict_multirateSample(X,D,sampFactors,N,f)

AD0 = reSample(N,D,sampFactors);
K = size(sampFactors,1);
if nargin < 5
    f = plotDictUsage(AD0,K,X);
else  
    f = plotDictUsage(AD0,K,X,f);
end

