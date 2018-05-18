%% Analyze texture
% Setup coordinates
fbasename = 'fista_fit_%i_%i.mat';
fname = sprintf(fbasename,0,0);
load(['D:\CHESS_data\al7075_311_space_ev_norm2\',fname])
step = 1;
num_bins = P.num_theta/step;

% Load Ring Fit
texture_az = zeros(num_bins,37,5,5);
texture_rad = zeros(num_bins,37,5,5);
signal = zeros(num_bins,37,5,5);

for load_num = 1:5
    for img_num = 0:184

        fbasename = 'fista_fit_%i_%i.mat';
        % Load Ring Fit
        fname = sprintf(fbasename,load_num-1,img_num);

        load(['D:\CHESS_data\al7075_311_space_ev_norm2\',fname])
        row = floor(img_num/5)+1;
        col = mod(img_num,5)+1;

        % Compute azimuthal variance ,radial variance, density at every theta
        for i = 1:num_bins
            ind = i:i+step-1;
            x_i = x_hat(:,ind,:,:);
            total = sum(x_i(:));

            % Azimuthal Variance
            for j = 1:numel(P.var_theta)
                x_ij = x_hat(:,ind,j,:);
                texture_az(i,row,col,load_num) = texture_az(i,row,col,load_num) + sum(x_ij(:))*(P.var_theta(j));
            end
            texture_az(i,row,col,load_num) = texture_az(i,row,col,load_num)/total;


            % Radial Variance

            for k = 1:numel(P.var_rad)
                x_ik = x_hat(:,ind,:,k);
                texture_rad(i,row,col,load_num) = texture_rad(i,row,col,load_num) + sum(x_ik(:))*(P.var_rad(k));
            end
            texture_rad(i,row,col,load_num) = texture_az(i,row,col,load_num)/total;


            % Total signal
            signal(i,row,col,load_num) = total;
        end
    end
end

disp('Texture')

%% Rebin signals
step = 128;
num_bins = P.num_theta/step;
texture_az_rebin = zeros(num_bins,37,5,5);
texture_rad_rebin = zeros(num_bins,37,5,5);
signal_rebin = zeros(num_bins,37,5,5);
for load_num = 1:5
    for img_num = 1:184
        row = floor(img_num/5)+1;
        col = mod(img_num,5)+1;
        for i = 1:num_bins
            ind = ((i-1)*step+1):i*step;
            for ii = ind
                texture_az_rebin(i,row,col,load_num) = texture_az_rebin(i,row,col,load_num) +...
                      texture_az(ii,row,col,load_num)*signal(ii,row,col,load_num)/total;
                texture_rad_rebin(i,row,col,load_num) = texture_rad_rebin(i,row,col,load_num) +...
                      texture_rad(ii,row,col,load_num)*signal(ii,row,col,load_num)/total;
                signal_rebin(i,row,col,load_num) = signal_rebin(i,row,col,load_num) +...
                                                   signal(ii,row,col,load_num);
            end
        end
    end
end

theta = 0:P.dtheta:(2*pi-P.dtheta);
theta_rebin = 0:P.dtheta*step:(2*pi-P.dtheta*step);
disp('Rebinning Complete')

%% Texture image
subplot(2,1,1)
imagesc(theta_rebin,1:5,squeeze(texture_az_rebin(:,36,1,:))')
colorbar()
subplot(2,1,2)
imagesc(theta,1:37,polar_image,[0 300])
colorbar()

%% Plot textures
figure
subplot(3,1,1)
bar(theta,texture_az)
title('\eta')

subplot(3,1,2)
bar(theta,texture_rad)
title('R')

subplot(3,1,3)
bar(theta,signal)
title('Signal')

