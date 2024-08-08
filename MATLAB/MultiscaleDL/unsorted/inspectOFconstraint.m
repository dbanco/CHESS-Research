function fig = inspectOFconstraint(u,v,Fx,Fy,Ft,window1,window2,t)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fig = figure
subplot(3,3,1)
imagesc(u(window1,window2,t))
title('u')
subplot(3,3,2)
imagesc(v(window1,window2,t))
title('v')


subplot(3,3,4)
imagesc(Fx(window1,window2,t))
title('Fx')
subplot(3,3,5)
imagesc(Fy(window1,window2,t))
title('Fy')

subplot(3,3,3)
imagesc(Ft(window1,window2,t)+Fy(window1,window2,t).*v(window1,window2,t)+Fx(window1,window2,t).*u(window1,window2,t))
title('Fxu + Fyv + Ft')

subplot(3,3,6)
imagesc(Fy(window1,window2,t).*v(window1,window2,t)+Fx(window1,window2,t).*u(window1,window2,t))
title('Fxu + Fyv')

subplot(3,3,7)
imagesc(Fx(window1,window2,t).*u(window1,window2,t))
title('Fx*u')

subplot(3,3,8)
imagesc(Fy(window1,window2,t).*v(window1,window2,t))
title('Fy*v')

subplot(3,3,9)
imagesc(Ft(window1,window2,t))
title('Ft')

end