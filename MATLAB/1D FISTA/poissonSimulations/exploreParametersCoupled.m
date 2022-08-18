function [mse,l1_norm,tv_penalty,awmv,fits,B] = exploreParametersCoupled(X_coupled,P,B)

[N,~,MM,T] = size(X_coupled);
A0 = unshifted_basis_vector_ft_stack_zpad(P);

mse = zeros(MM,1);
l1_norm = zeros(MM,1);
tv_penalty = zeros(MM,1);
awmv = zeros(T,MM);
fits = zeros(N,T,M);
for i = 1:MM
    for time = 1:T
        x = X_coupled(:,:,i,time);        
        fit = Ax_ft_1D(A0,x);
        fits(:,time,i) = fit;
        mse(i) = mse(i) + norm( B(:,time)-fit ).^2/norm(B(:,time))^2;
        l1_norm(i) = P.params.lambda1(time)*l1_norm(i) + sum(abs(x(:)));
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv(time,i) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end
    tv_penalty(i) = sum(abs(DiffPhiX_1D(X_coupled(:,:,i,:))),'all');
end
    
end

