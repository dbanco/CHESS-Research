function f = plotDictUsage(AD,K,X,f)
[~,~,KU] = size(AD);

if nargin < 4
    f = figure;
else
    f = figure(f);
end

for i = 1:KU
    subplot(KU/K,K,i)
    plot(AD(:,:,i)) 
%     if sum(squeeze(X(:,:,i,:)),'all')/sum(X(:)) > 0
%         title(sprintf('Usage: %1.2f',sum(squeeze(X(:,:,i,:)),'all')/sum(X(:)) ))
%     end
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
    p = 1+(f.Number-1)*400;
    f.Position =[p 1 400 1000];
end

end

