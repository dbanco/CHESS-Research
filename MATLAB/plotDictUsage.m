function f = plotDictUsage(D,rows,cols,X,vert,f)
[~,~,K] = size(D);

if nargin < 6
    f = figure;
else
    f = figure(f);
end
if nargin == 1
    rows = ceil(sqrt(K));
    cols = rows;
end
for i = 1:K
    if vert
        subplot(rows*cols,1,i)
    elseif cols >1
        subplot(cols,rows,i)
    else
        subplot(rows,cols,i)
    end
    plot(D(:,:,i)) 
    if sum(squeeze(X(:,:,i,:)),'all')/sum(X(:)) > 0
        title(sprintf('Usage: %1.2f',sum(squeeze(X(:,:,i,:)),'all')/sum(X(:)) ))
    end
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
    p = 1+(f.Number-1)*400;
    if cols < rows
        f.Position =[p 1 400 1000];
    else
        f.Position =[p 1 1000 400];
    end
end
end
