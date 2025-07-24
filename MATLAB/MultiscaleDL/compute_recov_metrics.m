function [error,rel_error] = compute_recov_metrics(outputs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    N = outputs.N;
    M = outputs.M;
    T = outputs.T;
    K = outputs.K;

    y = outputs.y;

    D = outputs.Dmin;
    X = outputs.Ymin;

    scales = outputs.scales;
    center = (M+1)/2;

    % Compute recons
    AD = reSampleCustomArrayCenter3(N,D,scales,center);
    AD = padarray(AD,[0 M-1 0],0,'post');
    ADf = fft2(AD);
    Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
    Yhat = gather(Yhat);
    % Compute error
    error = sum((squeeze(y)-Yhat).^2,'all');
    rel_error = error./sum(squeeze(y).^2,'all');
end