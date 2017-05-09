%% Setup
dr = 30;
radius = 370;
 
num_theta= 2048;
num_rad = 2*dr;

num_var_t = 15;
num_var_r = 10;

dtheta = 2*pi/num_theta;
drad = 1;

var_theta = linspace(dtheta,(pi/32),num_var_t).^2;
var_rad   = linspace(drad,3,num_var_r).^2;

A0_stack = unshifted_basis_matrix_stack(var_theta,...
                                           var_rad,...
                                           dtheta,...
                                           drad,...
                                           num_theta,... 
                                           num_rad);
                                       
A0_sum = squeeze(sum(sum(A0_stack,1),2));

result_path = fullfile('E:','CHESS_results','matlab');
backtrack_path = fullfile('E:','CHESS_results','backtrack');
error_dir = fullfile('E:','CHESS_results','matlab_error');
image_path = fullfile('E:','CHESS_data','matlab_polar_images');
%% Load results from files
% var_signal = zeros(num_var_t,num_var_r,5,41,5);
% rel_error = zeros(5,41,5);
load('spread_data1_5.mat')
% load_steps = 5;
% for step = 1:2
%     for img_num = 1:205
%         disp(['Load: ' num2str(step) '   Image: ' num2str(img_num)])
%         idx1 = floor((img_num-1)/5)+1;
%         idx2 = mod(img_num-1,5)+1;
%         try
%             file_name = ['backtrack_load_' num2str(step-1) '_img_' num2str(img_num-1) '.mat'];
%             file_path = fullfile(backtrack_path,file_name);
%             ringModel = load(file_path);
%             rel_error(step,idx1,idx2) = ringModel.rel_fit_error_new;
%             var_signal(:,:,step,idx1,idx2) = squeeze(sum(sum(ringModel.xhat_new,1),2)).*A0_sum;
% 
%         catch
%             file_name = ['ringModel_coefs_load_' num2str(step-1) '_img_' num2str(img_num-1) '.mat'];
%             file_path = fullfile(result_path,file_name);
%             ringModel = load(file_path);
%             var_signal(:,:,step,idx1,idx2) = squeeze(sum(sum(ringModel.xhat,1),2)).*A0_sum;
%             err_struc = load(fullfile(error_dir,... 
%             ['ringModel_error_load_',...
%             num2str(step-1), '_img_',...
%             num2str(img_num-1), '.mat']));
%             rel_error(step,idx1,idx2) = err_struc.rel_fit_error;
%         end
%     end
% end
% save('spread_data1_5.mat','var_signal','rel_error')
        
%% Plot results
cutoff_t = 3;
cutoff_r = 4;
load_steps = 5;
total_var = squeeze(sum(sum(var_signal(:,:,1:load_steps,:,:),1),2));
for i = 1:load_steps
    figure(1)
    high_var_theta = squeeze(sum(sum(var_signal(cutoff_t:end,:,1:load_steps,:,:),1),2))./total_var;
    all = high_var_theta(1:load_steps,3:40,:);
    mu = mean(all(:));
    sig = std(all(:));
    subplot(1,5, i)  
    imshow(squeeze(high_var_theta(i,:,:)),'DisplayRange',[0 1],'Colormap',jet)
    if( i == 2 )
        title('Theta Spread')
    end
    figure(2)
    high_var_rad = squeeze(sum(sum(var_signal(:,cutoff_r:end,1:load_steps,:,:),1),2))./total_var;
    all = high_var_rad(1:load_steps,3:40,:);
    mu = mean(all(:));
    sig = std(all(:));
    subplot(1,5, i)  
    imshow(squeeze(high_var_rad(i,:,:)),'DisplayRange',[0 1],'Colormap',jet)
    if(i ==2 )
        title('Radial Spread')
    end
   

    figure(3)
    all = rel_error(1:load_steps,3:40,:);
    mu = mean(all(:));
    sig = std(all(:));
    subplot(1,5, i)  
    imshow(squeeze(rel_error(i,:,:)),'DisplayRange',[0 0.5],'Colormap',jet)
    if( i ==2 )
        title('Fit Error')
    end
end
%% Load polar and fit image 
step = 1;
img_num = 30;

A0ft_stack = unshifted_basis_matrix_ft_stack(var_theta,...
                                           var_rad,...
                                           dtheta,...
                                           drad,...
                                           num_theta,... 
                                           num_rad);

disp(['Load: ' num2str(step-1) '   Image: ' num2str(img_num)])
file_name = ['polar_image_al7075_load_' num2str(step-1) '_img_' num2str(img_num) '.mat'];
file_path = fullfile(image_path,file_name);
load(file_path);
try
    file_name = ['backtrack_load_' num2str(step-1) '_img_' num2str(img_num) '.mat'];
    file_path = fullfile(backtrack_path,file_name);
    ringModel = load(file_path);
    fit = Ax_ft_2D(A0ft_stack,ringModel.xhat_new);
catch
    disp('No backtrack corrected file found')
    file_name = ['ringModel_coefs_load_' num2str(step-1) '_img_' num2str(img_num) '.mat'];
    file_path = fullfile(result_path,file_name);
    ringModel = load(file_path);
    fit = Ax_ft_2D(A0ft_stack,ringModel.xhat);
end

figure(4)
subplot(3,1,1) 
imshow(polar_image,'DisplayRange',[0 200],'Colormap',jet)
subplot(3,1,2)
imshow(fit,'DisplayRange',[0 200],'Colormap',jet)
subplot(3,1,3)
imshow(500*A0_stack(:,:,14,9),'DisplayRange',[0 200],'Colormap',jet)

%% View basis functions for individual spot

spot_org = polar_image(20:40,1204:1262);
spot_fit = fit(20:40,1204:1262);

% Count coefs (53780) a lot
dims = size(ringModel.xhat_new);
spot_coefs = ringModel.xhat_new(20:40,1204:1262,:,:);
ind0 = spot_coefs > 0;
l0_norm = sum(ind0(:));



%% Var coefs 
px = 6;
py = 3;

coef_sum = squeeze(sum(squeeze(sum(spot_coefs,1)),1));
coef_sum_w = coef_sum.*A0_sum;
% Show spot original
figure(5)
subplot(px,py,1) 
imshow(spot_org,'DisplayRange',[0 200],'Colormap',jet)
title(['coefs sum',' --- ','weighted coefs sum' ]);
% Show spot fit
figure(5)
subplot(px,py,2) 
imshow(spot_fit,'DisplayRange',[0 200],'Colormap',jet)


% Subplot coef values (fixed radial variance)
for i = 1:15
    subplot(px,py,i+2) 
    coef_slice = sum(spot_coefs(:,:,i,:),4);
    imshow(coef_slice,'DisplayRange',[0 1],'Colormap',jet)
    title([num2str(sum(coef_sum(i,:))),' --- ',num2str(sum(coef_sum_w(i,:))) ]);
end
% Subplot basis function (fixed azimuthal variance)
figure(9)
for i = 1:15
    subplot(px,py,i+2) 
    basis = A0_stack(:,:,i,10);
    basis_shift = shift2D(basis,10,29);
    imshow(basis_shift(1:21,1:59),'DisplayRange',[0 1],'Colormap',jet)
end

%% Rad coefs
px = 6;
py = 2;
coef_sum = squeeze(sum(squeeze(sum(spot_coefs,1)),1));
coef_sum_w = coef_sum.*A0_sum;



% Show spot original
figure(7)
subplot(px,py,1) 
imshow(spot_org,'DisplayRange',[0 200],'Colormap',jet)
title(['coefs sum',' --- ','weighted coefs sum' ]);
% Show spot fit
figure(7)
subplot(px,py,2) 
imshow(spot_fit,'DisplayRange',[0 200],'Colormap',jet)

% Subplot coef values (fixed azimuthal variance)
for i = 1:10
    subplot(px,py,i+2) 
    coef_slice = sum(spot_coefs(:,:,:,i),3);
    imshow(coef_slice,'DisplayRange',[0 1],'Colormap',jet)
    title([num2str(sum(coef_sum(:,i))),' --- ',num2str(sum(coef_sum_w(:,i))) ]);
end

% Subplot basis function (fixed azimuthal variance)
figure(8)
for i = 1:10
    subplot(px,py,i+2) 
    basis = A0_stack(:,:,7,i);
    basis_shift = shift2D(basis,10,29);
    imshow(basis_shift(1:21,1:59),'DisplayRange',[0 1],'Colormap',jet)
end
