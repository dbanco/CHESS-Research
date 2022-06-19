function f = checkDict_upDwnSample(X,D,U,f)
AD0 = upDwnSample(D,U);
if nargin < 4
    f = plotDictUsage(AD0,2*U+1,1,X,1);
else  
    f = plotDictUsage(AD0,2*U+1,1,X,1,f);
end

