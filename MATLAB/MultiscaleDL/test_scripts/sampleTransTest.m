
%% Test 1
close all

N = 91;
M = 51;
d = [zeros(1,15), 1:10,11,10:-1:1, zeros(1,15)];
centerM = (M+1)/2;
centerN = (N+1)/2;

c1 = 3;
c2 = 7;
% Resampling
du1 = uSample(d,c1);
dd1 = dSampleCenter(du1,c2,c1*centerM);
[~,Nud2,~] = size(dd1);
dd1_pad = padarray(dd1,[0,N-Nud2],0,'post');
dd1_padshift = circshift(dd1_pad, centerM-floor(centerM*c1/c2)-1);
% Transposed follow-through
ds2t = circshift(dd1_padshift,-(centerM-floor(centerM*c1/c2)-1));
du2t = uSampleCenter(ds2t,c2,c1*centerM);
dd2t = dSampleCenter(du2t,c1,c1*centerM);
dd_crop = dd2t(1:M);

% Passing sample X to transpose
x = zeros(N,1);
x(20) = 1;
x_shift = circshift(x,-(centerM-floor(centerM*c1/c2)-1));
x_up = uSampleCenter(x_shift,c2,c1*centerM);
x_down = dSampleCenter(x_up,c1,c1*centerM);
x_crop = x_down(1:M);

[N1,N2,Utotal] = size(dd1_padshift);
Nmin = min(N2,size(dd2t,2));

figure
subplot(3,4,1)
max1 = find(d==max(d));
plot(d,'o-')
title(['Initial at ',num2str(max1)])

subplot(3,4,2)
max2 = find(du1==max(du1));
plot(du1,'x-')
title(['Upsample at ',num2str(max2)])

subplot(3,4,3)
max3 = find(dd1==max(dd1));
plot(dd1,'x-')
title(['Downsample at ',num2str(max3)])

subplot(3,4,4)
max4 = find(dd1_padshift==max(dd1_padshift));
plot(dd1_padshift,'x-')
title(['Pad-shift at ',num2str(max4)])

subplot(3,4,5)
max5 = find(ds2t==max(ds2t));
plot(ds2t,'o-')
title(['Pad-shift back ',num2str(max5)])

subplot(3,4,6)
max6 = find(du2t==max(du2t));
plot(du2t,'o-')
title(['Upsample at ',num2str(max6)])

subplot(3,4,7)
max7 = find(dd2t==max(dd2t));
plot(dd2t,'o-')
title(['Downsample at ',num2str(max7)])

subplot(3,4,8)
max8 = find(dd_crop==max(dd_crop));
plot(dd_crop,'x-')
title(['Crop at ',num2str(max8)])

subplot(3,4,9)
max5 = find(x_crop==max(x_crop));
plot(x_crop,'o-')
title(['Pad-shift back ',num2str(max5)])

subplot(3,4,10)
max6 = find(x_up==max(x_up));
plot(x_up,'o-')
title(['Upsample at ',num2str(max6)])

subplot(3,4,11)
max7 = find(x_down==max(x_down));
plot(x_down,'o-')
title(['Downsample at ',num2str(max7)])

subplot(3,4,12)
max8 = find(x_crop==max(x_crop));
plot(x_crop,'x-')
title(['Crop at ',num2str(max8)])



% TEST CASE
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);
KJ = K*J;

d_test = rand(M,K);

[Ad,NormVals] = reSampleCustomArrayCenter(N,d_test,scales,centerM);

y_test = rand(1,N,KJ).*Ad;

Aty = reSampleTransCustomArrayCenter(M,y_test,scales,centerM,NormVals);

check1 = y_test(:)'*Ad(:)
check2 = d_test'*Aty(:)


% Adpad = padarray(Ad,[0,M-1, 0, 0],0,'post');


% c1 = 2;
% c2 = 3;
% % Resampling
% du3 = uSample(d,c1);
% dd3 = dSampleCenter(du3,c2,c1*center);
% dd3 = circshift(dd3,center-round(center*c1/c2));
% 
% % Transposed resampling
% du4 = uSampleCenter(d,c2,c1*center);
% dd4 = dSampleCenter(du3,c1,center);

% figure
% hold on
% plot(du1,'o-')
% plot(du2,'x')
% plot(du3,'*-')
% plot(du4,'s-')
% 
% figure
% hold on
% plot(dd1,'o-')
% plot(dd2,'x')
% plot(dd3,'*-')
% plot(dd4,'s-')
% 
% norm(dd1-dd4)
% norm(dd2-dd3)
