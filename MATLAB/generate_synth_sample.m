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

d_weight = max_az_var/5;

% Setup sample positions and 
positions = zeros(25,2);
dist = zeros(25,1);
center = [3,3];
k = 1;
for i = 1:5
    for j = 1:5
        positions(k,1) = i;
        positions(k,2) = j;
        dist(k) = sqrt((i-center(1))^2 + (j-center(2))^2);
        k = k + 1;
    end
end

% Produce 5x5 array of ring images
synth_sample = cell(5,5);
evar_image = zeros(5,5);
for i = 1:size(positions,1)
    row = positions(i,1);
    col = positions(i,2);
    
    % Generate diffraction spot parameters
    var_theta = rand(num_spots,1).^2*max_az_var + min_az_var + dist(i)*d_weight;
    var_rad = rand(num_spots,1).^2*max_rad_var + min_rad_var;
    theta_mean = (rand(num_spots,1)-0.5)*max_az_mean + min_az_mean;
    rad_mean = rand(num_spots,1)*max_rad_mean + min_rad_mean;
    amplitudes = rand(num_spots,1)*max_amp + min_amp;
    
    % Generate image
    B = zeros(num_rad,num_theta);
    for j = 1:num_spots
        B = B + amplitudes(j)*...
            gaussian_basis_wrap_2D( num_theta,dtheta,round(theta_mean(j)),var_theta(j),...
                                    num_rad,drad,round(rad_mean(j)),var_rad(j) );
    end
    % Add noise
    B = B + randn(num_rad,num_theta)*mean(amplitudes)/10;
    
    % Plot image
    % figure(2)
    % imshow(B,'DisplayRange',[0 300],'Colormap',jet)
    
    % Save sample info to structure
    sample.image = B;
    sample.var_theta = var_theta;
    sample.var_rad = var_rad;
    sample.theta_mean = theta_mean;
    sample.rad_mean = rad_mean;
    sample.amplitudes = amplitudes;
    % Peak function generating parameters
    sample.P.max_amp = max_amp;
    sample.P.min_amp = min_amp;
    sample.P.max_az_mean = max_az_mean;
    sample.P.min_az_mean = min_az_mean;
    sample.P.max_az_var = max_az_var;
    sample.P.min_az_var = min_az_var;
    sample.P.max_rad_mean = max_rad_mean;
    sample.P.min_rad_mean = min_rad_mean;
    sample.P.max_rad_var = max_rad_var;
    sample.P.min_rad_var = min_rad_var;

    %Compute expected variance
    a = sum(amplitudes(:));
    b = amplitudes.*var_theta;
    expected_var = sum(b(:))/a;
    sample.expected_var = expected_var;
    evar_image(row,col) = expected_var;
    
    % Add sample to 5x5 cell array 
    synth_sample{row,col} = sample;
    im_name = sprintf('F:\\CHESS_data\\synth_data\\polar_image_0_%i.mat',i);
    polar_image = B;
    
    % Save out image files
    save(im_name,'polar_image')
end

%save('F:\CHESS_data\synth_data\synth_data.mat','synth_sample')

%% View Expected Variance
figure(2)
image(evar_image,'CDataMapping','scaled')
colorbar()