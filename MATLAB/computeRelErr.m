function relErr = computeRelErr(B,X_hat,A0)

[N,T] = size(B);
relErr = zeros(T,1);
for t = 1:T
    b = B(:,t);
    x = squeeze(X_hat(:,:,t));
    fit = Ax_ft_1D(A0,x);
%     figure(22)
%     hold on
%     plot(b)
%     plot(fit)
%     pause()
    relErr(t) = norm(fit-b)/norm(b);
end
