%% Compute center of mass
top_dir = 'E:\MMPAD_data';
dset_name = 'ring%i_zero';
prefix = 'mmpad_img_%i.mat';

R = 4;
N = 546;
COM = zeros(R,N);
for ring = 1:R
    for img = 1:N
        ring_dir = sprintf(dset_name,ring);
        img_file = sprintf(prefix,img);
        load(fullfile(top_dir,ring_dir,img_file))
        
        [n,m] = size(polar_image);
        radial_img = sum(polar_image,2);
        radial_coord = (1:n)';
        total_mass = sum(radial_img);
        COM(ring,img) = sum(radial_img.*radial_coord./total_mass);
    end
end


%% Load AWMV
datadir = 'C:\Users\dpqb1\Desktop\present_figures';
fileBase = 'ring%i_zero_couple_fit_data.mat';
outName =  'ring%i_CoM_AWMV.png';
for ring = 1:4
    fig_out = figure(ring);
    load(fullfile(datadir,sprintf(fileBase,ring)))
    
    % Rescale center of Mass
    minAWMV = min(awmv_az(select_ind,:));
    maxAWMV = max(awmv_az(select_ind,:));
    scaleCOM = (COM(ring,:)-min(COM(ring,:)))/(max(COM(ring,:))-min(COM(ring,:)));
    scaleCOM = (scaleCOM + minAWMV/maxAWMV)*maxAWMV;
    hold on
    title(sprintf('Ring %i',ring))
    plot(awmv_az(select_ind,:),'LineWidth',1.5)
    plot(scaleCOM,'LineWidth',1.5)
    legend('AWMV','scaled CoM')
    saveas(fig_out,fullfile(datadir,sprintf(outName,ring)))
end

%% PLot all awmv
datadir = 'C:\Users\dpqb1\Desktop\present_figures';
fileBase = 'ring%i_zero_couple_fit_data.mat';
outName =  'allAWMV.png';
fig_out = figure(5);
for ring = 1:4
    
    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(awmv_az(select_ind,:),'LineWidth',1.5)
    ylabel('AWMV')
    xlabel('time')

end
legend('ring 1','ring 2','ring 3','ring 4')
saveas(fig_out,fullfile(datadir,outName))

%% Plot all CoM
outName =  'allCoM.png';
fig_out = figure(6);
for ring = 1:4
    
    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(COM(ring,:),'LineWidth',1.5)
    ylabel('CoM')
    xlabel('time')

end
legend('ring 1','ring 2','ring 3','ring 4')
saveas(fig_out,fullfile(datadir,outName))