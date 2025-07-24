function sigma = estimate_noise(y,window)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

half_win = (window-1)/2;

y = squeeze(y);
[N,T] = size(y);

variances = zeros(N-window+1,T);

mean_filt = ones(window,1)/window;

means_window = convn(y,mean_filt,'valid');


for t = 1:T
    i = 1;
    for n = 1+half_win:N-half_win
        n1 = n-half_win;
        n2 = n+half_win;
        variances(i,t) = mean((y(n1:n2)-means_window(i)).^2);
        i = i + 1;
    end
end

sigma = sqrt(mode(variances(:)));
