% Test stationary wavelet transform with triangle scaling function
yn = gaus_linear_osc_signal(0.005);

phi = [0 0.25 0.5 0.25];
psi = [0 -0.25 0.5 -0.25];

[swa,swd] = swt(yn(1,:,1),3,phi,psi);
swc = swt(yn(1,:,1),3,phi,psi);

plot(swc')
figure
plot(swd')

