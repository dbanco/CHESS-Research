
%% Test 1
close all
d = [zeros(1,15), 1:10,11,10:-1:1, zeros(1,15),zeros(1,60)];
center = 26;

%22*2/3 = 14.667 -> 15

c1 = 4;
c2 = 5;
% Resampling
du1 = uSample(d,c1);
dd1 = dSampleCenter(du1,c2,c1*center);
ds1 = circshift(dd1,center-round(center*c1/c2));
% Transposed follow-through
ds2t = circshift(ds1,round(center*c1/c2)-center);
du2t = uSampleCenter(ds2t,c2,c1*center);
dd2t = dSampleCenter(du2t,c1,c1*center);
% Transposed resampling on d
% ds2 = circshift(d,round(center*c1/c2)-center);
% du2 = uSampleCenter(ds2,c2,c1*center);
% dd2 = dSampleCenter(du2,c1,center);

figure
plot(d,'o-')

figure
subplot(2,3,1)
plot(du1,'x-')
subplot(2,3,2)
plot(dd1,'x-')
subplot(2,3,3)
plot(ds1,'x-')

subplot(2,3,4)
plot(ds2t,'o-')
subplot(2,3,5)
plot(du2t,'o-')
subplot(2,3,6)
hold on
plot(dd2t,'o-')
% plot(d,'s-')
plot(ds1,'x-')




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
