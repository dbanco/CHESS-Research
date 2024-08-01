function [Yhat1,Yhat2] = fig3decomp(outputs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = outputs.N;
M = outputs.M;
D = gather(outputs.D);
X = gather(outputs.X);
y = gather(outputs.y);
scales = outputs.scales;
center = (M+1)/2;

% generateFiguresToy1zpad_center(figDir,outputs,suffix,[4,4]);

D = outputs.D;
X = outputs.X;
scales = outputs.scales;
N = outputs.N;
M = outputs.M;
y = outputs.y;
K = outputs.K;
center = (M+1)/2;
Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
Utotal = sum(Uarray);

AD = reSampleCustomArrayCenter(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);



Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));
Yhat = unpad(Yhat,M-1,'pre');

X1 = X;
X1(:,:,9:16,:) = 0;
Yhat1 = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X1)),3),'symmetric'));
Yhat1 = unpad(Yhat1,M-1,'pre');

X2 = X;
X2(:,:,1:8,:) = 0;
Yhat2 = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X2)),3),'symmetric'));
Yhat2 = unpad(Yhat2,M-1,'pre');
end