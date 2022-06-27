function checkDict(X,D,U,f)
K = size(D,3);
AD0 = decimate(D,U);
if nargin < 4
    f = plotDictUsage(AD0,K,X);
else  
    f = plotDictUsage(AD0,K,X,f);
end

end

