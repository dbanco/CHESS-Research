function f = plotDictUsage(AD,K,X,f)
[~,~,KU] = size(AD);

if nargin < 4
    f = figure;
else
    f = figure(f);
end

for i = 1:KU
    subplot(2,ceil(KU/2),i)
    plot(AD(:,:,i),'Linewidth',1) 
%     set(gca,'visible','off');
%     ylim([0 1])
%     if sum(squeeze(X(:,:,i,:)),'all')/sum(X(:)) > 0
%         title(sprintf('Usage: %1.2f',sum(squeeze(X(:,:,i,:)),'all')/sum(X(:)) ))
%     end
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
%     p = 1+(f.Number)*400;
    f.Position(3:4) =[1800 400];
end

end

