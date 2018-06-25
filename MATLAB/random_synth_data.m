%% Randomized spot example
num_spots = 200;
% Ring sampling parameters
ring_width = 20;
num_theta= 2048;
num_omega = 360;
num_rad = 2*ring_width+1;
dtheta = 2*pi/num_theta;
drad = 1;
domega = 2*pi/num_omega;

% Amplitude
min_amp = 100;
max_amp = 300;

%Azimuthal mean and variance
min_az_mean = 1;
max_az_mean = num_theta-1;
min_az_var = dtheta^2;
max_az_var = (pi/150)^2;

%Omega mean and variance
min_om_mean = 1;
max_om_mean = num_omega-1;
min_om_var = domega^2;
max_om_var = (pi/150)^2;

%Radial mean and variance
min_rad_mean = ring_width;
max_rad_mean = 0.7;
min_rad_var = drad^2;
max_rad_var = 4;

var_theta = rand(num_spots,1).^2*max_az_var + min_az_var;
var_rad = rand(num_spots,1).^2*max_rad_var + min_rad_var;
var_omega = rand(num_spots,1).^2*max_om_var + min_om_var;
theta_mean = (rand(num_spots,1)-0.5)*max_az_mean + min_az_mean;
rad_mean = rand(num_spots,1)*max_rad_mean + min_rad_mean;
omega_mean = rand(num_spots,1)*max_om_mean + min_om_mean;
amplitudes = rand(num_spots,1)*max_amp + min_amp;

B = zeros(num_rad,num_theta,num_omega);
for i = 1:num_spots
    fprintf('Spot %i \n',i)
    B = B + amplitudes(i)*gaussian_basis_wrap_3D( num_theta,dtheta,round(theta_mean(i)),var_theta(i),...
                                                  num_rad,drad,round(rad_mean(i)),var_rad(i),...
                                                  num_omega,domega,round(omega_mean(i)),var_omega(i));   

end
B = B+randn(num_rad,num_theta,num_omega)*mean(amplitudes)/10;

%%
figure(2)
for i = 1:num_omega
    imshow(squeeze(B(:,:,i)),'DisplayRange',[0 300],'Colormap',jet)
    pause(0.05)
end

