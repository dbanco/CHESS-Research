saveDir = 'D:\CHESS_data\simulated_two_spot_1D\';
mkdir(saveDir)

%% Randomized spot example

% Basis function variance parameters
P.num_var_t = 300;
P.dtheta = 1;

% Ring sampling parameters
num_theta= 401;
dtheta = 1;

az_mean1 = 201 - (0:0.5:45);
az_mean2 = 201 + (0:0.5:45);
az_std_both = 10;
amplitudes = [1,1];

num_ims = numel(az_mean2);

synth_sample = cell(num_ims,1);
VDF = cell(num_ims,1);
evar_az = zeros(num_ims,1);
allB = [];

figure(222)
[ha, pos] = tight_subplot(11,10,[.005 .005],[.01 .01],[.01 .01]); 
k = 1;
for i = 1:num_ims

%         VDF{i} = histogram(az_std_both(j),P.var_theta');

        % Generate image
        B = zeros(num_theta,1);
        B = B + gaussian_basis_1D_norm2( num_theta, az_mean1(i), az_std_both^2);
                                   
        B = B + gaussian_basis_1D_norm2( num_theta, az_mean2(i), az_std_both^2);
                                          
        % Add noise
%         B = B + randn(num_rad,num_theta)/100;

        % Plot image
        axes(ha(i));

        plot(B)
        
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
% 
%         % Add sample to (n1*n2)x1 cell array 
%         synth_sample{i} = sample;
        im_name = [saveDir,sprintf('polar_vector_%i.mat',i)];
        polar_vector = B;

        % Save out image files
        save(im_name,'polar_vector')
end

% save([saveDir,'synth_data.mat'],'synth_sample','VDF')

%% View Expected Variance
% figure(2)
% image(evar_image,'CDataMapping','scaled')
% colorbar()