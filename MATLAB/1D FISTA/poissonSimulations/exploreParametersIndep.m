function [mse,l1_norm,awmv,fits,B] = exploreParametersIndep(X_indep,P,B)

[N,~,M,T] = size(X_indep);
A0 = unshifted_basis_vector_ft_stack_zpad(P);

mse = zeros(T,M);
l1_norm = zeros(T,M);
awmv = zeros(T,M);
fits = zeros(N,T,M);
for i = 1:M
    for time = 1:T
        x = X_indep(:,:,i,time);
        fit = Ax_ft_1D(A0,x);
        fits(:,time,i) = fit;
        % Objectives
        l1_norm(time,i) = sum(abs(x(:)));
        mse(time,i) = norm(fit-B(:,time));
        % AWMV
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv(time,i) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end
end
    
end

