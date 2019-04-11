saveDir = 'D:\CHESS_data\simulated_data_nooverlap_60\';
mkdir(saveDir)
%% Randomized spot example
num_spots = 60;
% Ring sampling parameters
ring_width = 20;
num_theta= 2048;
num_rad = 2*ring_width+1;
dtheta = 1;
drad = 1;
% Amplitude
min_amp = 100;
max_amp = 300;

%Azimuthal mean and variance
min_az_mean = 1;
max_az_mean = 2047;
start_center_az_std = 2;
end_center_az_std = 25;
width_az_std = 1;

%Radial mean and variance
center_rad_mean = 20;
width_rad_mean = 5;
start_center_rad_std = 1;
end_center_rad_std = 4;
width_rad_std = 0.3;

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.dtheta = 1;
P.drad = 1;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t);
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r);

% Setup sample positions and 
positions = zeros(100,1);
dist = zeros(100,1);
for i = 1:100
        positions(i) = i;
        dist(i) = sqrt( (i-1)^2 );
end

% Produce 5x5 array of ring images
synth_sample = cell(100,1);
VDF = cell(100,1);
evar_az = zeros(100,1);
evar_rad = zeros(100,1);

theta_means = (rand(num_spots,1)-0.5)*max_az_mean + min_az_mean;
rad_means = rand(num_spots,1)*width_rad_mean + center_rad_mean;
amplitudes = rand(num_spots,1)*max_amp + min_amp;

for i = 1:size(positions,1)
% for i = [1,50,100]
    % Generate diffraction spot parameters
    std_theta = randn(num_spots,1).*width_az_std + start_center_az_std*(1-dist(i)/100) + end_center_az_std*(dist(i)/100);
    std_rad = randn(num_spots,1).*width_rad_std + start_center_rad_std*(1-dist(i)/100) + end_center_rad_std*(dist(i)/100); 

    VDF{i} = hist3([std_theta,std_rad],'Ctrs',{P.var_theta',P.var_rad'});
    figure(1)
    imagesc(VDF{i})
    
    % Generate image
    B = zeros(num_rad,num_theta);
    for j = 1:num_spots
        B = B + amplitudes(j)*...
            gaussian_basis_wrap_2D_norm2( num_theta,dtheta,round(theta_means(j)),std_theta(j)^2,...
                                          num_rad,drad,round(rad_means(j)),std_rad(j) );
    end
    % Add noise
    B = B + randn(num_rad,num_theta)*mean(amplitudes)/100;
    
    % Plot image
    figure(2)
    uplim = max(B(:));
    imshow(B,'DisplayRange',[0 uplim],'Colormap',jet)
    
    % Save sample info to structure
    sample.image = B;
    sample.std_theta = std_theta;
    sample.std_rad = std_rad;
    sample.theta_mean = theta_means;
    sample.rad_mean = rad_means;
    sample.amplitudes = amplitudes;
    sample.vdf = VDF{i};
    
    %Compute expected variance
    total = sum(amplitudes(:));
    az_signal = amplitudes.*std_theta;
    rad_signal = amplitudes.*std_rad;
    expected_var_az = sum(az_signal(:))/total;
    expected_var_rad = sum(rad_signal(:))/total;
    evar_az(i) = expected_var_az;
    evar_rad(i) = expected_var_rad;
    
    % Add sample to 100x1 cell array 
    synth_sample{i} = sample;
    im_name = [saveDir,sprintf('polar_image_%i.mat',i)];
    polar_image = B;
    
    % Save out image files
    save(im_name,'polar_image')
    
    % Perturb spots
    theta_mean_deltas = randn(num_spots,1)*2;
    rad_mean_deltas = randn(num_spots,1)*0.25;
    amplitude_deltas = randn(num_spots,1)*8;
    
    theta_means = theta_means + theta_mean_deltas;
    rad_means = rad_means + rad_mean_deltas;
    amplitudes = amplitudes + amplitude_deltas;
    pause()
end
save([saveDir,'synth_data.mat'],'synth_sample','VDF')

%% View Expected Variance
% figure(2)
% image(evar_image,'CDataMapping','scaled')
% colorbar()