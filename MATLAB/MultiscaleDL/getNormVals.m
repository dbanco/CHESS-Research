function NormVals = getNormVals(Dout)
[~,~,KU] = size(Dout);
NormVals = zeros(KU,1);
for i = 1:KU
    NormVals(i) = norm(Dout(:,:,i));
end