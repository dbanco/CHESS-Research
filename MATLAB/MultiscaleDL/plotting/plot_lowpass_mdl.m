%% View lowpass filter for up/down sampling
close all
m = 5;
M = 2*m-1;

f = zeros(1,M);
f(1:m) = linspace(1/m,1,m);
f(m+1:end) = f(m-1:-1:1);
f = f/norm(f);

n = 3;
N = 2*n-1;

f2 = zeros(1,N);
f2(1:n) = linspace(1/n,1,n);
f2(n+1:end) = f2(n-1:-1:1);
f2 = f2/norm(f2);

f3  = conv(f,f2);
NM = numel(f3);

topDir = 'C:\Users\dpqb1\Desktop\dissertation_figures';

fig1 = figure;
plot(padarray(f',(NM-M)/2),'-or','Linewidth',2)
ylim([0 1])
set(gca,'XColor','none','YColor','none')
saveas(fig1,fullfile(topDir,'phi1.png'))

fig2 = figure;
plot(padarray(f2',(NM-N)/2),'-og','Linewidth',2)
ylim([0 1])
set(gca,'XColor','none','YColor','none')
saveas(fig2,fullfile(topDir,'phi2.png'))

fig3 = figure;
plot(f3,'-ob','Linewidth',2)
ylim([0 1])
set(gca,'XColor','none','YColor','none')
saveas(fig3,fullfile(topDir,'phi1phi2.png'))

% set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])

% str1 = sprintf('\\phi_%i',m);
% str2 = sprintf('\\phi_%i',n);
% str3 = sprintf('\\phi_%i \\ast \\phi_%i',m,n);
% legend(str1,str2,str3,'FontSize',14)
