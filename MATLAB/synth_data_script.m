%% Powder example
% Ring sampling parameters
ring_width = 20;
num_theta= 2048;
num_rad = 2*ring_width+1;
dtheta = 2*pi/num_theta;
drad = 1;

% Basis function variance parameters
var_theta = (pi)^2;
var_rad   = (4)^2;
theta_mean = round(num_theta/2);
rad_mean = ring_width;

B = gaussian_basis_wrap_2D( num_theta,dtheta,theta_mean,var_theta,...
                                 num_rad,drad,rad_mean,var_rad );
B = B+B.*randn(41,2048)/20;     
figure(1)
imshow(B,'DisplayRange',[0 1],'Colormap',jet)

%% Spot example

% Ring sampling parameters
ring_width = 20;
num_theta= 2048;
num_rad = 2*ring_width+1;
dtheta = 2*pi/num_theta;
drad = 1;

% Basis function variance parameters
var_theta = [pi/300, pi/320, pi/350].^2;
var_rad   = [4 3.2 2.2].^2;
theta_mean = [round(num_theta/2-47), round(num_theta/2), round(num_theta/2+35)];
rad_mean = [ring_width,ring_width-1,ring_width+2];

B = zeros(num_rad,num_theta);
for i = 1:3
    B = B + gaussian_basis_wrap_2D( num_theta,dtheta,theta_mean(i),var_theta(i),...
                                     num_rad,drad,rad_mean(i),var_rad(i) );
    B = B+B.*randn(41,2048)/20;   
end
figure(2)
imshow(B,'DisplayRange',[0 1],'Colormap',jet)