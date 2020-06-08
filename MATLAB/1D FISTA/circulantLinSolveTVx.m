function [x_hat] = circulantLinSolveTVx( A0ft_stack,b,yk,vk,z1k,z2k,u1k,u2k,x_n,params )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg
%

bnormsq = sum(b(:).^2);
b_ft = fft(b);

rho = params.rho*bnormsq;
rho2 = params.rho2*bnormsq;

tm1 = numel(x_n{1}) > 1;
tp1 = numel(x_n{2}) > 1;

if tm1 && tp1
    const_ft = fft(rho1*(yk-vk)+rho2*(x_n{1}+z1k-u1k)+rho2*(x_n{2}-z2k+u2k));
elseif tm1
    const_ft = fft(rho1*(yk-vk)+rho2*(x_n{1}+z1k-u1k));
elseif tp1
    const_ft = fft(rho1*(yk-vk)+rho2*(x_n{2}-z2k+u2k));
end
    
r_ft = (A0ft_stack'.').*repmat(b_ft,[1,20])/bnormsq + const_ft;
      
coef = diag(A0ft_stack*r_ft.')./((rho + (tm1+tp1)*rho2) +...
       diag(A0ft_stack*A0ft_stack'));
x_ft = (r_ft - repmat(coef,[1,20]).*(A0ft_stack'.'))/...
       (rho + (tm1+tp1)*rho2);
x_hat = real(ifft(x_ft));
end

