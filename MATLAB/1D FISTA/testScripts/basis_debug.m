%%
a = gaussian_basis_wrap_1D_norm2(10,0,1)
b = gaussian_basis_wrap_1D_norm2(10,5.5,1)
% c = gaussian_basis_wrap_1D_norm2(11,0,1)
figure(1) 
hold on
plot(a,'o-')
plot(b,'o-')
% plot(c)
%%
a = gaussian_basis_wrap_1D_norm2(11,0,1)
b = gaussian_basis_wrap_1D_norm2(11,6,1)
% c = gaussian_basis_wrap_1D_norm2(11,0,1)
figure(2) 
hold on
plot(a,'o-')
plot(b,'o-')
% plot(c)