function [output] = goldenScaleSearch(y,N,D0,X_true,U,denLim)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
gr = (sqrt(5) + 1)/2;
x = [1, 4-3/gr, 1+3/gr, 4];
probes = [1 0 0 4;1 0 0 1];
for i = 2:3
    [probes(1,i),probes(2,i)] = fareyApprox(x(i),denLim);
end
err = zeros(4,1);
for i = 1:4
    err(i) = evalProbe(y,probes(:,i),N,D0,X_true,U);
end
[minErr,ind] = min(err);
output = probes(:,ind);

while minErr > 1e-4
    if err(2) < err(3)
        x(4) = x(3);
        probes(:,4) = probes(:,3);
    else
        x(1) = x(2);
        probes(:,1) = probes(:,2);
    end
    x(2) = x(4)-(x(4)-x(1))/gr;
    x(3) = x(1)+(x(4)-x(1))/gr;
    for i = 2:3
        [probes(1,i),probes(2,i)] = fareyApprox(x(i),denLim);
    end
    for i = 1:4
        err(i) = evalProbe(y,probes(:,i),N,D0,X_true,U);
    end
    [minErr,ind] = min(err);
    output = probes(:,ind);

    figure(3)
    hold on
    plot(probes(1,:)./probes(2,:),err,'-o')
    probes
end

end

