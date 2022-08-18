function [x_hat, dcx_hat] = circulantLinSolve_DC( A0ft_stack,b,yk,vk,params,dcyk,dcvk)
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg

% bnormsq = sum(b(:).^2);
N = numel(b);
b_ft = fft(b);
c_ft = fft(ones(N,1));
K = size(A0ft_stack,2);
rho = params.rho1;

r_ft = (A0ft_stack'.').*repmat(b_ft-c_ft*dcyk,[1,K]) + rho*fft(yk-vk);
coef = diag(A0ft_stack*r_ft.')./...
      (rho + diag(A0ft_stack*A0ft_stack'));
x_ft = (r_ft - repmat(coef,[1,K]).*(A0ft_stack'.') )/rho;
x_hat = real(ifft(x_ft));

N = numel(b);
Axt = Ax_ft_1D(A0ft_stack,yk)';
dcx_hat = (b'*ones(N,1) + Axt*ones(N,1) + rho*dcvk)/(N + rho);

% if params.plotProgress 
%     figure(2)
%     hold off
%     fit = Ax_ft_1D(A0ft_stack,x_hat) + dcx_hat;
%     plot(b)
%     hold on
%     plot(fit)
%     legend('b','fit')
% end

end

