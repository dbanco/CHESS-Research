%% Randomized spot example
num_spots = 100;
% Ring sampling parameters
ring_width = 20;
num_theta= 2048;
num_rad = 2*ring_width+1;
dtheta = 2*pi/2048;
drad = 1;

% Amplitude
min_amp = 100;
max_amp = 300;
%Azimuthal mean and variance
min_az_mean = 1;
max_az_mean = 2047;
min_az_var = dtheta^2;
max_az_var = (pi/150)^2;

%Radial mean and variance
min_rad_mean = ring_width;
max_rad_mean = 0.7;
min_rad_var = drad^2;
max_rad_var = 4;

var_theta = rand(num_spots,1).^2*max_az_var + min_az_var;
var_rad = rand(num_spots,1).^2*max_rad_var + min_rad_var;
theta_mean = (rand(num_spots,1)-0.5)*max_az_mean + min_az_mean;
rad_mean = rand(num_spots,1)*max_rad_mean + min_rad_mean;
amplitudes = rand(num_spots,1)*max_amp + min_amp;

B = zeros(num_rad,num_theta);
for i = 1:num_spots
    B = B + amplitudes(i)*gaussian_basis_wrap_2D( num_theta,dtheta,round(theta_mean(i)),var_theta(i),...
                                     num_rad,drad,round(rad_mean(i)),var_rad(i) );   

end
B = B+randn(num_rad,num_theta)*mean(amplitudes)/10;

figure(2)
imshow(B,'DisplayRange',[0 300],'Colormap',jet)

