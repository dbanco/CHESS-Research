function [D, Y, optinf, obj, relErr,output,minObj] = cbpdndlScaleSearch(D0,S,lambda,U,denLim,opt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
probes = [1, 2, 3, 4;
          1, 1, 1, 1];
Jiters = 200;
opt.MaxMainIter = Jiters;
opt.Verbose = 0;
states = struct([]);
% Begin solutions for the 4 probe points
for i = 1:4
    fprintf('Updating Probe %i: ',i)
    [D, Y, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, S, lambda, opt,probes(1,i),probes(2,i),U);
    fprintf('relErr = %0.3f\n', relErr)
    states(i).D = D;
    states(i).Y = Y;
    states(i).opt = optinf;
    states(i).obj = obj;    
end

[minObj,ind] = min([states.obj]);
% opt.Y0 = states(ind).Y;
% D0 = states(ind).D;
output = probes(:,ind);

while relErr > 0.01
    if states(2).obj < states(3).obj
        probes(:,4) = probes(:,3);
        states(4) = states(3);
        len1 = probes(1,4)/probes(2,4) - probes(1,2)/probes(2,2);
        len2 = probes(1,2)/probes(2,2) - probes(1,1)/probes(2,1);
        if len2 < len1
            probes(:,3) = probes(:,2) + probes(:,4);
            fprintf('Updating Probe %i: ',3)
            [D, Y, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, S, lambda, opt,probes(1,3),probes(2,3),U);
            fprintf('relErr = %0.3f\n', relErr)
            states(3).D = D;
            states(3).Y = Y;
            states(3).opt = optinf;
            states(3).obj = obj;
            if states(3).obj  < minObj
                ind = 3;
                minObj = states(3).obj ;
                output = probes(:,3);
            end
            if probes(2,3) > denLim; break; end
        else
            probes(:,3) = probes(:,2);
            states(3).obj  = states(2).obj ;
            probes(:,2) = probes(:,1) + probes(:,2);
            fprintf('Updating Probe %i: ',2)
            [D, Y, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, S, lambda, opt,probes(1,2),probes(2,2),U);
            fprintf('relErr = %0.3f\n', relErr)
            states(2).D = D;
            states(2).Y = Y;
            states(2).opt = optinf;
            states(2).obj = obj;
            if states(2).obj < minObj
                ind = 2;
                minObj = states(2).obj;
                output = probes(:,2);
            end
            if probes(2,2) > denLim; break; end
        end
    else
        probes(:,1) = probes(:,2);
        states(1)= states(2);
        len1 = probes(1,4)/probes(2,4) - probes(1,3)/probes(2,3);
        len2 = probes(1,3)/probes(2,3) - probes(1,1)/probes(2,1);
        if len1 < len2
            probes(:,2) = probes(:,1) + probes(:,3);
            fprintf('Updating Probe %i as %i/%i ',2,probes(1,2),probes(2,2))
            [D, Y, optinf, obj,relErr] = cbpdndl_cg_multirate(D0, S, lambda, opt,probes(1,2),probes(2,2),U);
            fprintf('relErr = %0.3f, obj = %2.3f\n', relErr, obj)
            states(2).D = D;
            states(2).Y = Y;
            states(2).opt = optinf;
            states(2).obj = obj;
            if states(2).obj < minObj
                ind = 2;
                minObj = states(2).obj;
                output = probes(:,2);
            end
            if probes(2,2) > denLim; break; end
        else
            probes(:,2) = probes(:,3);
            states(2) = states(3);
            probes(:,3) = probes(:,3) + probes(:,4);
            fprintf('Updating Probe %i: ',3)
            [D, Y, optinf, obj] = cbpdndl_cg_multirate(D0, S, lambda, opt,probes(1,3),probes(2,3),U);
            fprintf('relErr = %0.3f\n', relErr)
            states(3).D = D;
            states(3).Y = Y;
            states(3).opt = optinf;
            states(3).obj = obj;
            if states(3).obj < minObj
                ind = 3;
                minObj = states(3).obj;
                output = probes(:,3);
            end
            if probes(2,3) > denLim; break; end
        end
    end 
    figure(3)
    hold on
    plot(probes(1,:)./probes(2,:),[states.obj],'-o')
    probes
end

fprintf('Final output is %i/%i\n',output(1),output(2))    
D = states(ind).D;
Y = states(ind).Y;
optinf = states(ind).opt;
obj = states(ind).obj;

    
end

