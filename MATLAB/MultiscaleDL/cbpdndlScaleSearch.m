function [D, Y, optinf, obj, relErr,output,minObj,prbCount] = cbpdndlScaleSearch(D0,S,lambda,U,denLim,opt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
probes = [1, 2, 3, 4;
          1, 1, 1, 1];
opt.Verbose = 0;
states = struct([]);

% Begin solutions for the 4 probe points
fprintf('Probes: % 2i % 2i % 2i % 2i \n',probes(1,:))
fprintf('        % 2i % 2i % 2i % 2i \n',probes(2,:))
for i = 1:4
    fprintf('Updating Probe %i as %i/%i ',i,probes(1,i),probes(2,i))
    [D, Y, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, S, lambda, opt,probes(1,i),probes(2,i),U);
    fprintf('relErr = %0.3f,  obj = %0.3f\n', relErr, obj)
    states(i).D = D;
    states(i).Y = Y;
    states(i).opt = optinf;
    states(i).obj = obj;    
end

objPlot = figure;
stem(probes(1,:)./probes(2,:),[states.obj],'o')
hold on


[minObj,ind] = min([states.obj]);
% opt.Y0 = states(ind).Y;
% D0 = states(ind).D;
output = probes(:,ind);
prbCount = 4;
while relErr > 0.01
    prbCount = prbCount + 1;
    if states(2).obj < states(3).obj
        probes(:,4) = probes(:,3);
        states(4) = states(3);
        len1 = probes(1,4)/probes(2,4) - probes(1,2)/probes(2,2);
        len2 = probes(1,2)/probes(2,2) - probes(1,1)/probes(2,1);
        if len2 < len1
            probes(:,3) = probes(:,2) + probes(:,4);
            [states,ind,minObj,output] = updateState(3,D0, S, lambda, opt,probes,states,U,ind,minObj,output,objPlot);
            if probes(2,3) > denLim; break; end
        else
            probes(:,3) = probes(:,2);
            states(3).obj  = states(2).obj ;
            probes(:,2) = probes(:,1) + probes(:,2);
            [states,ind,minObj,output] = updateState(2,D0, S, lambda, opt,probes,states,U,ind,minObj,output,objPlot);
            if probes(2,2) > denLim; break; end
        end
    else
        probes(:,1) = probes(:,2);
        states(1)= states(2);
        len1 = probes(1,4)/probes(2,4) - probes(1,3)/probes(2,3);
        len2 = probes(1,3)/probes(2,3) - probes(1,1)/probes(2,1);
        if len1 < len2
            probes(:,2) = probes(:,1) + probes(:,3);
            [states,ind,minObj,output] = updateState(2,D0, S, lambda, opt,probes,states,U,ind,minObj,output,objPlot);
            if probes(2,2) > denLim; break; end
        else
            probes(:,2) = probes(:,3);
            states(2) = states(3);
            probes(:,3) = probes(:,3) + probes(:,4);
            [states,ind,minObj,output] = updateState(3,D0, S, lambda, opt,probes,states,U,ind,minObj,output,objPlot);
            if probes(2,3) > denLim; break; end
        end
    end 
end

fprintf('Final output is %i/%i\n',output(1),output(2))
D = states(ind).D;
Y = states(ind).Y;
optinf = states(ind).opt;
obj = states(ind).obj;

end

function [states,ind,minObj,output] = updateState(i,D0, S, lambda, opt,probes,states,U,ind,minObj,output,objPlot);
fprintf('Probes: % 2i % 2i % 2i % 2i \n',probes(1,:))
fprintf('        % 2i % 2i % 2i % 2i \n',probes(2,:))
fprintf('Updating Probe %i as %i/%i ',i,probes(1,i),probes(2,i))
[D, Y, optinf, obj, relErr] = cbpdndl_cg_multirate(D0, S, lambda, opt,probes(1,i),probes(2,i),U);
fprintf('relErr = %0.3f,  obj = %0.3f\n', relErr, obj)
states(i).D = D;
states(i).Y = Y;
states(i).opt = optinf;
states(i).obj = obj;
if states(i).obj < minObj
    ind = i;
    minObj = states(i).obj;
    output = probes(:,i);
end
if nargin > 8
    objPlot = figure(objPlot);
    stem(probes(1,i)/probes(2,i),obj,'o')
    objPlot.Position = [1 750 500 300];
end
end
