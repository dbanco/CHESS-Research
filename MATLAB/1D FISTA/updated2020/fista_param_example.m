tic
res = [1.53 0.9 0.80 0.73 0.65 0.46 0.42 0.40 0.39 0.39 0.39 0.39];
L1_norm = [4 110 219 309 462 2521 7435 15941 35481 46944 46962 46962];


figure(1)
plot(res,L1_norm,'o-')
xlabel('Residual')
ylabel('||x||_1')
title('seismic tomography paper')

figure(2)
loglog((res),(L1_norm),'o-')
xlabel('Residual')
ylabel('||x||_1')
title('seismic tomography paper')
toc