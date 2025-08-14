%% Plot example multiscale dicitonary
figDir = "C:\Users\dpqb1\Documents\MCDL Paper";
N = 100;
M = 80;
center = floor((N+1)/2);
eta_dom = linspace(-3.5,3.5,M);
D = [zeros(1,10),sinc(eta_dom) + 0.25];
D = D./norm(D);
scales = {};
scales{1} = genRationals([0;1],[1;1],4,10, 1/4);
%scales{2} = genRationals([0;1],[1;1],4,10, 1/4);

AD = reSampleCustomArrayCenter3(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);

f12j = cell(6,1);
f13j = cell(6,1);
for i = 0:5
    f12j{i+1} = figure(120+i);
    plot(AD(:,:,6-i),'Linewidth',2)
    xlim([0 N])
    ylim([0 0.5])
    f12j{i+1}.Position = [150  500-50*i  200 100];
    set(gca,'FontSize',12)
    % xt = [0,50,80,100];
    % xtlabels = string(xt);
    % xticks(xt)
    saveas(f12j{i+1},fullfile(figDir,['atom_K1_J',num2str(i+1),'.png']))
end