%% Load MMPAD sequence

for r = 1
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
lb = -1;
ub = 1;
n = 10;
[psi,~] = gauswavf(lb,ub,n,2);
figure(11)
plot(psi)

bank = dwtfilterbank('SignalLength',M,...
                   'Wavelet','Custom',...
                   'FilterType','Analysis',...
                   'CustomWaveletFilter',psi',...
                   'CustomScalingFilter',psi');
     
figure(2)
waveletScale = zeros(T,1);
for t = 1:541
    x = X1d(:,t);  
    Xout = filterbankAna(x,bank,5);
    waveletScale(t) = computeWaveletScale(Xout);
end

% Delta waveletScale
waveletScale = waveletScale - waveletScale(1);
figure(2)
hold on
plot(waveletScale)

end
ylabel('\Delta waveletScale')
xlabel('time')


