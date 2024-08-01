%% Load MMPAD sequence
close all
waveType = {'bior4.4','bior3.9','bior1.1','haar'};
rings = {'\{004\}','\{021\}','\{112\}','\{020\}'};

r = 1
data_dir = ['D:\MMPAD_data_nr1\ring', num2str(r), '_zero'];
T = 546;

for t = 1:T
	load(fullfile(data_dir,['mmpad_img_',num2str(t),'.mat']))
    if t == 1
        [N,M] = size(polar_image);
        X = zeros(N,M,T);
    end
    X(:,:,t) = polar_image;
end

X1d = squeeze(sum(X,1));

% 2-layer Filter bank decomp on mmpad

bank = dwtfilterbank('SignalLength',M,...
                   'Wavelet',waveType{1},...
                   'FilterType','Analysis');
               
figure(2)
waveletScale = zeros(T,1);
for t = 100
    x = X1d(:,t);  
    nn = 4;
    Xout = filterbankAna(x,bank,nn);
    waveletScale(t) = computeWaveletScale(Xout);
    Xcomp = filterbankVis(x,Xout,bank);
    
    figure(44)
    for ii = 1:2^nn
        subplot(4,4,ii)
        plot(Xcomp(:,ii))
    end
end


figure(55)
plot(x)
%%
% Delta waveletScale
waveletScale = waveletScale - waveletScale(1);
figure(2)
hold on
plot(waveletScale)


ylabel('\Delta waveletScale')
xlabel('time')
c_leg = legend(rings,'Location','Best');
ctitle = get(c_leg,'Title');
set(ctitle,'String','Ring indices \{hkl\}')

