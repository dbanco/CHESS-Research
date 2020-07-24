function [x_hat] = circulantLinSolve_2D( A0ft_stack,b,yk,vk,params )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg

[N1,N2,K1,K2] = size(A0ft_stack);
N = N1*N2;
K = K1*K2;
A0 = zeros(N,K);
x_ft_array = zeros(size(A0ft_stack));
bnormsq = sum(b(:).^2);
b1_ft = reshape(fft2(b),[N,1]);

rho = params.rho1*bnormsq;

term2 = rho*fft2(yk-vk);
term2NK = zeros(N,K);

k = 1;
for k1 = 1:K1
    for k2 = 1:K2
        A0(:,k) = reshape(A0ft_stack(:,:,k1,k2),[N,1]);
        term2NK(:,k) = reshape(term2(:,:,k1,k2),[N,1]);
        k = k + 1;
    end
end

r_ft = conj(A0).*repmat(b1_ft,[1,K]) + term2NK;
coef = sum(A0.*r_ft,2)./ (rho + sum(A0.*A0,2));
x_ft = (r_ft - repmat(coef,[1,K]).*conj(A0) )/rho;

k = 1;
for k1 = 1:K1
    for k2 = 1:K2
        x_ft_array(:,:,k1,k2) = reshape(x_ft(:,k),[N1,N2]);
        k = k + 1;
    end
end

x_hat = real(ifft2(x_ft_array));
end

