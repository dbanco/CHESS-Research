function [error,rel_error,of_penalty, hs_penalty, D_error] = compute_recov_metrics(outputs,dataset,sigma,lambda2,lambda3,HSiters)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    N = outputs.N;
    M = outputs.M;
    T = outputs.T;
    K = outputs.K;
    J = outputs.J;

    y = outputs.y;

    D = outputs.D;
    Y = outputs.Y;

    scales = outputs.scales;
    center = (M+1)/2;

    % Compute recons
    AD = reSampleCustomArrayCenter3(N,D,scales,center);
    AD = padarray(AD,[0 M-1 0],0,'post');
    ADf = fft2(AD);
    Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(Y)),3),'symmetric')),M-1,'pre');
    Yhat = gather(Yhat);
    % Compute error
    error = sum((squeeze(y)-Yhat).^2,'all');
    rel_error = error./sum(squeeze(y).^2,'all');

    % Compute OF, HS penalties
    [Uvel,Vvel,Fx,Fy,Ft] = computeHornSchunkDictPaperLS(Y,K,outputs.Uvel,outputs.Vvel,lambda3/lambda2,HSiters);
    [of_penalty, hs_penalty] = HSobjectivePaper(Fx,Fy,Ft,Uvel,Vvel,K,lambda3/lambda2);
    
    % Dict Error
    [~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigma,dataset);
    if ~isscalar(Dtrue)
        
        [D_perm,~] = align_third_dim_and_shift(D, Dtrue);

        figure
        subplot(2,1,1)
        plot(D_perm(1,:,1))
        hold on
        plot(Dtrue(1,:,1))
        subplot(2,1,2)
        plot(D_perm(1,:,2))
        hold on
        plot(Dtrue(1,:,2))

 
        % Compute errors on recovered X and D 
        D_error = sqrt(sum((D_perm-Dtrue).^2,'all'))/sqrt(sum((Dtrue).^2,'all'));
        
    end


end