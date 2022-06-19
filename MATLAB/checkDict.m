function checkDict(X,D,U)
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
AD0 = decimate(D,U);
plotDictUsage(AD0,5,1,X,1);

end

