saveDir = 'D:\CHESS_data\simulated_two_spot_1D_anomaly';

%% Randomized spot example
num_ims = 20;
num_spots = 2;

% Basis function variance parameters
P.dtheta = 1;

% Ring sampling parameters
num_theta= 180;
dtheta = 1;

synth_sample = cell(num_ims,1);
VDF = cell(num_ims,1);
evar_az = zeros(num_ims,1);
allB = [];

theta_means = [45,125];

amplitudes = [1;1];

% theta_stds1 = [2:16];
% theta_stds2 = [8:20];

t = linspace(0,1,num_ims);
b = 10*(1-t) + 5*t;
amp1 = (1-t)*2 + t*0.1;
amp2 = (1-t)*3 + t*3;


theta_stds1 = [8*ones(1,num_ims/2),12*ones(1,num_ims/2)];
theta_stds2 = [8*ones(1,num_ims/2),12*ones(1,num_ims/2)];

figure(222)
[ha, pos] = tight_subplot(15,10,[.005 .005],[.01 .01],[.01 .01]); 
kk = 1;

noise_std = [0:0.03:0.30];
for k = 1:numel(noise_std)
    for i = 1:num_ims

    %         VDF{i} = histogram(az_std_both(j),P.var_theta');

        % Generate image
        B = zeros(num_theta,1);
        Bn = zeros(num_theta,1);
        B = B + amplitudes(1)*gaussian_basis_1D( num_theta, theta_means(1)+randn(1)/2, theta_stds1(i)^2);
        B = B + amplitudes(2)*gaussian_basis_1D( num_theta, theta_means(2)+randn(1)/2, theta_stds2(i)^2);

        % Add noise
        for j = 1:numel(B)
            norm(B);
            norm(randn(1)*noise_std(k));
            Bn(j) = B(j) + randn(1)*noise_std(k);
        end

        % Don't remove negative component
%         Bn = max(Bn, 0);

%             B = B + randn(num_rad,num_theta)/100;
%             rel_error = norm(Bn-B)/norm(B)
%             snr = norm(B)^2/norm(Bn-B)^2

        % Plot image
        if k == 1
            try
                axes(ha(kk));
                kk = kk + 1;
                plot(Bn)
            catch

            end
        end

        % Save sample info to structure
        sample.image = B;
        sample.std_theta = [theta_stds1(i); theta_stds2(i)];
        sample.theta_mean = [theta_means(1); theta_means(2)];
        sample.amplitudes = [1;1];
        sample.vdf = VDF{i};

        %Compute expected variance
        total = sum(amplitudes(:));
        az_signal = amplitudes.*sample.std_theta;
        expected_var_az = sum(az_signal(:))/total;
        evar_az(i) = expected_var_az;
        
        % Add sample to (n1*n2)x1 cell array 
        synth_sample{i} = sample;
        im_name = [saveDir,sprintf('_%i/',k),sprintf('polar_vector_%i.mat',i)];
        mkdir([saveDir,sprintf('_%i/',k)])
        polar_vector = Bn;
% 
%             % Save out image files
            save(im_name,'polar_vector')
    end
end
figure(1)
plot(evar_az)

save([[saveDir,sprintf('_%i',k)],'\synth_data.mat'],'synth_sample','VDF')

%% View Expected Variance
% figure(2)
% image(evar_image,'CDataMapping','scaled')
% colorbar()