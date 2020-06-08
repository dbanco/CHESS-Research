function [x_hat] = circulantLinSolveTVx( A0ft_stack,b,yk,vk,zk1,zk2,u1k,u2k,params,neighbors )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg
%

rho = params.rho;
rho2 = params.rho2;
if params.time == 1
    factor = -1;
elseif neighbors == 1
    factor = 1;
else
    factor = 0;
end

b_ft = fft(b);
y_ft = fft(yk);
v_ft = fft(vk);
z1_ft = fft(zk1);
z2_ft = fft(zk2);
u1k_ft = fft(u1k);
u2k_ft = fft(u2k);

r_ft = (A0ft_stack'.').*repmat(b_ft,[1,20]) + rho*(y_ft-v_ft) -...
     rho2*(z1_ft-u1k_ft) + rho2*(z2_ft-u2k_ft);

coef = diag(A0ft_stack*r_ft.')./((rho + factor*rho2) + diag(A0ft_stack*A0ft_stack'));
x_ft = (r_ft - repmat(coef,[1,20]).*(A0ft_stack'.'))/(rho + factor*rho2);
x_hat = real(ifft(x_ft));
end

