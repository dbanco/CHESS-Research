imagesc(img(:,:,1))
hold on
center = [1944,3560];
r = 1215;
theta = linspace(0,2*pi,300);

x = r*cos(theta) + center(1);
y = r*sin(theta) + center(2);

plot(x,y)
