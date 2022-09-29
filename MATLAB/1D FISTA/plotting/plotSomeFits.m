function figOut = plotSomeFits(X_hat,B,A0,Time)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
figOut = figure;
i = 1;
TT = numel(Time);
for t = Time
    b = B(:,t);
    x = squeeze(X_hat(:,:,t));
    fit = Ax_ft_1D(A0,x); 

    subplot(ceil(sqrt(TT)),ceil(sqrt(TT)),i)
    hold on
    plot(b)
    plot(fit)
    i = i + 1;
    title(num2str(norm(b-fit)/norm(b)))
end

