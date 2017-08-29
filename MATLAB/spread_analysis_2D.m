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
                                       
result_path = fullfile('D:','CHESS_results','fista_fit_results');
image_path = fullfile('D:','CHESS_data','matlab_polar_images');
%% Get results
%  var_signal = zeros(num_var_t,num_var_r,5,41,5);
%  rel_error = zeros(5,41,5);
load('spread_data_load1.mat')

load_steps = 5;
for step = 2:load_steps
    for img_num = 1:205
        disp(['Load: ' num2str(step) '   Image: ' num2str(img_num)])
        idx1 = floor((img_num-1)/5)+1;
        idx2 = mod(img_num-1,5)+1;
        
        file_name = ['fista_out_load_' num2str(step-1) '_img_' num2str(img_num-1) '.mat'];
        file_path = fullfile(result_path,file_name);
        load(file_path);
        rel_error(step,idx1,idx2) = err(end);
        var_signal(:,:,step,idx1,idx2) = squeeze(sum(sum(x_hat,1),2));
    end
end
save('spread_data_all.mat','var_signal','rel_error')
        
%% Plot results
cutoff_t = 2;
cutoff_r = 3;
load_steps = 5;
total_var = squeeze(sum(sum(var_signal(:,:,1:load_steps,:,:),1),2));
for i = 1:load_steps
    figure(1)
    high_var_theta = squeeze(sum(sum(var_signal(cutoff_t:end,:,1:load_steps,:,:),1),2))./total_var;
    subplot(1,5, i)  
    imshow(squeeze(high_var_theta(i,:,:)),'DisplayRange',[0 1],'Colormap',jet)
    if( i == 2 )
        title('Theta Spread')
    end
    
    figure(2)
    high_var_rad = squeeze(sum(sum(var_signal(:,cutoff_r:end,1:load_steps,:,:),1),2))./total_var;
    subplot(1,5, i)  
    imshow(squeeze(high_var_rad(i,:,:)),'DisplayRange',[0 0.5],'Colormap',jet)
    if(i ==2 )
        title('Radial Spread')
    end
   
    figure(3)
    subplot(1,5, i)  
    imshow(squeeze(rel_error(i,:,:)),'DisplayRange',[0 1],'Colormap',jet)
    if( i ==2 )
        title('Fit Error')
    end
end