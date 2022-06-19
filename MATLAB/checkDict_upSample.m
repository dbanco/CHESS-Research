function checkDict_upSample(X,D,U)
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
AD0 = upSample(D,U);
plotDictUsage(AD0,U,1,X,1);

end

