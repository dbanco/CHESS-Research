% Create figure showing the spreading of diffraction signal
close all
clear all
%% Load data
topDir = 'C:\Users\dpqb1\Documents\Data\ring1_zero';
y = loadMMPAD2D(1,0,topDir);
[N,M,T] = size(y);

%% Concatenate time steps and show

%Concatenate a bunch series of time-steps
t = [1,20,30,35,40,45,50,55,60,80,100,120];

subset = 135:200;
maxs = zeros(numel(t),1);
fullImg = zeros(numel(subset),N*numel(t));
dottedLine = zeros(numel(t),1);
for i = 1:numel(t)
    inds = (1:N) + N*(i-1);
    yt = y(:,:,t(i))';
    maxs(i) = max(yt(:));
    yt = yt/max(yt(:));
    fullImg(:,inds) = yt(subset,:);
end

imagesc(fullImg)
hold on
for i = 1:numel(t)-1
    plot([1 1]*N*i,...
         [1 numel(subset)],...
        'w--')
end

axis image
set(gca, 'Visible','off')