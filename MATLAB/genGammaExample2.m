% rng('default'); % For reproduceability.
% pdTrue = GeneralizedGamma(1.4, 1.0, 1.6);
% n = 10000;
% sample = pdTrue.drawSample(n);
% pdEstimated = GeneralizedGamma();
% pdEstimated.fitDist(sample)

% Test #1: Does the implementation reproduce results from the literature? 
r=1; % mode location
eta = linspace(-1,1,5);
beta=3/2+eta;
vartheta=1; % like variance, larger = heavier tail, lower mode



fig1 = figure('position', [100 100 600 300]);
legend_str = cell(numel(r),1);
hold on

for i = 1:numel(beta)
    c = r;
    m = beta(i);
    lambda = 1/vartheta;
    pdTrue = GeneralizedGamma(lambda,c,m);
    x = [0:0.01:8];
    f = pdTrue.pdf(x);
    F = pdTrue.cdf(x);
    ny = median(f);
    ni = find(f==median(f));
    nx = [x(ni)*1.1, x(ni)]/8;
    ny = [ny*1.1, ny];
    plot(x, f);
    n2 =sprintf( '%0.2f',beta(i));
    string_beta = ['\beta=',n2,''];
    legend_str{i} = string_beta;
%     annotation('textarrow',nx,ny,'String', string_r);
    
end
legend(legend_str,'location','best','Fontsize',12')
ylabel('\alpha','Fontsize',16');
xlabel('\theta','Fontsize',16');
ylim([0 0.6])

ha = annotation('textbox',[0.6 0.7 0.12 0.16], 'Interpreter', 'latex');
n1 =sprintf( '%0.2f',r);
n3 =sprintf( '%0.2f',vartheta);
sbox = ['$r=',n1,'$\\ $\vartheta=',n3,'$'];
% sbox = ['$r=',n1,'$\\ $\beta=',n2,'$\\ $\vartheta=',n3,'$'];
set(ha, 'String', sbox,'Fontsize',12');
box off