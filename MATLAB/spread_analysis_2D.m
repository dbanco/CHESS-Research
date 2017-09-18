load('spread_alphap_10_betap_1.mat')

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