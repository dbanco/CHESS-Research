function plotDataSeq(y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y = squeeze(y);
T = size(y,2);

figure
for t = 1:T
    plot(y(:,t))
    pause(0.2)
end



end

