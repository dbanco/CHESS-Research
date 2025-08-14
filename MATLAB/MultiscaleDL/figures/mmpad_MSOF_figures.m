% Show simulated data

load("C:\Users\dpqb1\Documents\Outputs\mmpad_ring1_optFlow_K2_M261\output_j_5_10_5_lam1_2.00e-02_lam2_5.00e-01_lam3_1.00e-03.mat")
suffix = '';
figDir = "C:\Users\dpqb1\Documents\Dissertation\mmpad_MSOF";

N = outputs.N;
M = outputs.M;
D = gather(outputs.D);  
X = gather(outputs.X);
y = gather(outputs.y);
scales = outputs.scales;
center = (M+1)/2;

%% Coefficients in time
figure(1)
X_eta = squeeze(sum(X,3));
imagesc(X_eta)

%% Coefficients in scale

X_sigma = squeeze(sum(X,2));

f2a = figure(2);
imagesc(X_sigma(1:8,:))
f2a.Position = [100  500  500 120];
ylabel('\sigma')
xlabel('t')

f2b = figure(3);
imagesc(X_sigma(9:16,:))
f2b.Position = [100  300  500 120];
ylabel('\sigma')
xlabel('t')

%%
% saveas(f2a,fullfile(figDir,['vdf1','.png']))
% saveas(f2b,fullfile(figDir,['vdf2','.png']))

%% Plot multiscale dict atoms
D = outputs.D;
N = outputs.N;
M = outputs.M;
y = outputs.y;
X = outputs.X;
scales = outputs.scales;
center = (M+1)/2;

AD = reSampleCustomArrayCenter(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
Yhat = gather(Yhat);

% f12 = figure(12);
% plot(AD(:,:,8),'Linewidth',2)
% xlim([0 261])
% f12.Position = [150  500  340 110];
% 
% f13 = figure(13);
% plot(AD(:,:,16),'Linewidth',2)
% xlim([0 261])
% f13.Position = [500  500  340 110];

f12j = cell(8,1);
f13j = cell(8,1);
for i = 0:7
    f12j{i+1} = figure(120+i);
    plot(AD(:,:,8-i),'Linewidth',2)
    xlim([0 261])
    ylim([0 1])
    f12j{i+1}.Position = [150  500-50*i  340 110];
    saveas(f12j{i+1},fullfile(figDir,['atom_K1_J',num2str(i+1),'.png']))

    f13j{i+1} = figure(130+i);
    plot(AD(:,:,16-i),'Linewidth',2)
    xlim([0 261])
    ylim([0 0.33])
    f13j{i+1}.Position = [500  500-50*i  340 110];
%     saveas(f13j{i+1},fullfile(figDir,['atom_K2_J',num2str(i+1),'.png']))
end

%%
close all
fig_pos = [100  100  300 200];
fig_pos2 = fig_pos+[400 0 0 0];
azimuth = 167;
elevation = 45;

ff1=figure;
waterfall(Yhat')
ff1.Position = fig_pos;
xlim([0,261])
ylim([0,120])
clim([0,0.9])
view(azimuth, elevation);
xticks([61 161,261]); 
yticks([0, 120]); 
xticklabels({'200','100','0'});
yticklabels({'0', '120'});


ff2=figure;
waterfall(squeeze(y)')
ff2.Position = fig_pos2;
% xlim([0,261])
% ylim([0,120])
clim([0,0.9])
view(azimuth, elevation);
xticks([61 161,261]); 
yticks([0, 120]); 
xticklabels({'200','100','0'});
yticklabels({'0', '120'});

ff3=figure;
waterfall(squeeze(y)')
ff3.Position = fig_pos2;
xlim([0,261])
ylim([0,120])
clim([0,0.9])
view(azimuth, elevation);
xticks([61 161,261]); 
yticks([0, 120]); 
xticklabels({'200','100','0'});
yticklabels({'0', '120'});
colorbar()

% Colorbar
colorbar();


saveas(ff1,'waterfall_recon.png')
saveas(ff2,'waterfall_data.png')
saveas(ff3,'waterfall_colorbar.png')


ff4=figure;
waterfall(log(squeeze(y)')+5)
ff4.Position = fig_pos2;
xlim([0,261])
ylim([0,120])
clim([log(0.01)+5,log(0.9)+5])
view(azimuth, elevation);
xticks([61 161,261]); 
yticks([0, 120]); 
xticklabels({'200','100','0'});
yticklabels({'0', '120'});
colorbar()

%% Surface plot
ff5 = figure;
y = loadMMPAD1D(3,'copy-shift');
Z = squeeze(y(:,110:210,:))';  % Size: [numY, numX]

% Generate matching X and Y coordinates
[numY, numX] = size(Z);
[X, Y] = meshgrid(1:numX, 1:numY);  % match X and Y to Z

% Plot as surface
surf(X, Y, Z, 'EdgeColor', 'none');  % No mesh lines for a cleaner surface

% Apply view settings
ff5.Position = fig_pos2;
% xlim([0, 261])
% ylim([0, 120])
clim([0, 0.9])
view(azimuth, elevation);

% Customize ticks
% xticks([61 161 261]);
% yticks([0 120]);
% xticklabels({'200', '100', '0'});
% yticklabels({'0', '120'});