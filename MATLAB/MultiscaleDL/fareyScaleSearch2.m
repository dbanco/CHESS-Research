function [output,reconErr] = fareyScaleSearch(y,N,D0,X_true,U,denLim)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
probes = [1, 2, 3, 4;
          1, 1, 1, 1];
for i = 1:4
    err(i) = evalProbe(y,probes(:,i),N,D0,X_true,U);
end
[reconErr,ind] = min(err);
output = probes(:,ind);
while reconErr > 1e-2
    if err(2) < err(3)
        probes(:,4) = probes(:,3);
        err(4) = err(3);
        len1 = probes(1,4)/probes(2,4) - probes(1,2)/probes(2,2);
        len2 = probes(1,2)/probes(2,2) - probes(1,1)/probes(2,1);
        if len2 < len1
            probes(:,3) = probes(:,2) + probes(:,4);
            err(3) = evalProbe(y,probes(:,3),N,D0,X_true,U);
            if err(3) < reconErr
                if probes(2,3) > denLim; return; end
                reconErr = err(3);
                output = probes(:,3);
            end
        else
            probes(:,3) = probes(:,2);
            err(3) = err(2);
            probes(:,2) = probes(:,1) + probes(:,2);
            err(2) = evalProbe(y,probes(:,2),N,D0,X_true,U);
            if err(2) < reconErr
                if probes(2,2) > denLim; return; end
                reconErr = err(2);
                output = probes(:,2);
            end
        end
    else
        probes(:,1) = probes(:,2);
        err(1) = err(2);
        len1 = probes(1,4)/probes(2,4) - probes(1,3)/probes(2,3);
        len2 = probes(1,3)/probes(2,3) - probes(1,1)/probes(2,1);
        if len1 < len2
            probes(:,2) = probes(:,1) + probes(:,3);
            err(2) = evalProbe(y,probes(:,2),N,D0,X_true,U);
            if err(2) < reconErr
                if probes(2,2) > denLim; return; end
                reconErr = err(2);
                output = probes(:,2);
            end
        else
            probes(:,2) = probes(:,3);
            err(2) = err(3);
            probes(:,3) = probes(:,3) + probes(:,4);
            err(3) = evalProbe(y,probes(:,3),N,D0,X_true,U);
            if err(3) < reconErr
                if probes(2,3) > denLim; return; end
                reconErr = err(3);
                output = probes(:,3);
            end
        end
    end 
    figure(3)
    hold on
    plot(probes(1,:)./probes(2,:),err,'-o')
end
end

