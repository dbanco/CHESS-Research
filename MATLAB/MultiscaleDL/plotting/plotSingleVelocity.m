function fig = plotSingleVelocity(X,K,opt,figNum)

if nargin < 4
    figNum = 11;
end

fig = figure(figNum);
[u,v,~,~,~]  = computeHornSchunkDict(X,K,opt.Smoothness,opt.HSiters);
% True solution
subplot(3,1,1)
imagesc(squeeze(sum(v,3))')
title('v summed in time')
ylabel('scale')
xlabel('shift')

subplot(3,1,2)
imagesc(squeeze(sum(u,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')

subplot(3,1,3)
imagesc(squeeze(sum(X,2)))
title('Atom usage in time')
ylabel('scale')
xlabel('time')
end