saveDir = 'D:\CHESS_data\simulated_one_spot\';
mkdir(saveDir)

%% Randomized spot example

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.dtheta = 1;
P.drad = 1;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t);
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r);

% Ring sampling parameters
ring_width = 15;
num_theta= 80;
num_rad = 2*ring_width+1;
dtheta = 1;
drad = 1;

az_mean = 40;
rad_mean = 15;
az_std = 1:0.5:15;
rad_std = 0.5:0.25:5;

amplitudes = [1];

% Produce 5x5 array of ring images
n_az = numel(az_std); 
n_rad = numel(rad_std); 

synth_sample = cell(n_rad*n_az,1);
evar_az = zeros(n_rad*n_az,1);
allB = [];

figure(222)
[ha, pos] = tight_subplot(n_rad,n_az,[.005 .005],[.01 .01],[.01 .01]); 

for i = 1:n_rad
    for j = 1:n_az
        figure(1)

        % Generate image
        B = zeros(num_rad,num_theta);
        B = B + gaussian_basis_2D( num_theta, az_mean,  az_std(j)^2,...
                                   num_rad,   rad_mean, rad_std(i)^2 );
 
        % Add noise
%         B = B + randn(num_rad,num_theta)/100;

        % Plot image
%         ii = ((i-1)*n_az + j);
%         axes(ha(ii));
%         uplim = max(B(:));
%         imshow(B,'DisplayRange',[0 uplim],'Colormap',jet)
        
%         % Save sample info to structure
%         sample.image = B;
%         sample.std_theta = [az_std_both(j); az_std_both(j)];
%         sample.std_rad = [rad_std; rad_std];
%         sample.theta_mean = [az_mean1; az_mean2(i)];
%         sample.rad_mean = [rad_mean; rad_mean];
%         sample.amplitudes = [1;1];
%         sample.vdf = VDF{i};
% 
%         %Compute expected variance
%         total = sum(amplitudes(:));
%         az_signal = amplitudes.*sample.std_theta;
%         rad_signal = amplitudes.*std_rad;
%         expected_var_az = sum(az_signal(:))/total;
%         expected_var_rad = sum(rad_signal(:))/total;
%         evar_az(i) = expected_var_az;
%         evar_rad(i) = expected_var_rad;

        % Add sample to (n1*n2)x1 cell array 
%         synth_sample{i} = sample;
        im_name = [saveDir,sprintf('polar_image_%i_%i.mat',i,j)];
        polar_image = B;

        % Save out image files
        save(im_name,'polar_image')
    end
end

%% View Expected Variance
% figure(2)
% image(evar_image,'CDataMapping','scaled')
% colorbar()