%% Powder example
% Ring sampling parameters
ring_width = 20;
num_theta= 200;
num_rad = 2*ring_width+1;
num_omega= 20;

dtheta = 2*pi/num_theta;
drad = 1;
domega = 2*pi/num_omega;
% Basis function variance parameters
var_theta = (pi/30)^2;
var_rad   = (4)^2;
var_omega = (pi/10)^2;

theta_mean = 0;
rad_mean = 0;
omega_mean = 0;

B = gaussian_basis_wrap_3D_norm2(num_theta,dtheta,theta_mean,var_theta,...
                                 num_rad,drad,rad_mean,var_rad,...
                                 num_omega,domega,omega_mean,var_omega);
B = B+B.*randn(num_rad,num_theta,num_omega)/20;

figure(1)
for i = 1:num_omega
    imshow(squeeze(B(:,:,i)),'DisplayRange',[0 0.1],'Colormap',jet)
    pause(0.1)
end
%% Spot example

% Ring sampling parameters
ring_width = 20;
num_theta= 500;
num_rad = 2*ring_width+1;
num_omega = 50;
dtheta = 2*pi/num_theta;
drad = 1;
domega = 2*pi/num_omega;

% Basis function variance parameters
var_theta = [pi/30, pi/32, pi/35].^2;
var_rad   = [4 3.2 2.2].^2;
var_omega = [pi/30, pi/32, pi/35].^2;
theta_mean = [round(num_theta/2-47), round(num_theta/2), round(num_theta/2+35)];
rad_mean = [ring_width,ring_width-1,ring_width+2];
omega_mean = [round(num_omega/2-47), round(num_omega/2), round(num_omega/2+35)];

B = zeros(num_rad,num_theta,num_omega);
for i = 1:3
    B = B + gaussian_basis_wrap_3D_norm2( num_theta,dtheta,theta_mean(i),var_theta(i),...
                                          num_rad,drad,rad_mean(i),var_rad(i),...
                                          num_omega,domega,omega_mean(i),var_omega(i));
    B = B+B.*randn(num_rad,num_theta,num_omega)/20;   
end
figure(2)
for i = 1:num_omega
    imshow(squeeze(B(:,:,i)),'DisplayRange',[0 0.1],'Colormap',jet)
    pause(0.1)
end