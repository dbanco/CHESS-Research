function f = plotDictionary(D,rows,cols)
D = squeeze(D);
[~,n2] = size(D);
f = figure;
if nargin == 1
    rows = ceil(sqrt(n2));
    cols = rows;
end
for i = 1:n2
    subplot(rows,cols,i)
    plot(D(:,i)) 
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
end
end

