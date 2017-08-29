function plot_circle( center, radius )
%plot_circle Summary of this function goes here
%   Detailed explanation goes here

theta = linspace(0,2*pi,5000);
x = radius*cos(theta) + center(1);
y = radius*sin(theta) + center(2);
plot(x,y,'w-')

end

