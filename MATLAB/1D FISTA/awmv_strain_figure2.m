% MMPAD metadata
mmpad_dims = [396,265];
rings = {'020','112','021','004'};
d_space = [1.27773,1.24845,1.23272,1.17138];
two_theta = [6.89132,7.05316,7.14328,7.5179];

pixel_side = 150e-6;
gap_width = 0.75e-3;

detec_dist = 4.6;
detect_angle = 14.6;

% Compute circumeference of rings
radius = 4.6*tan(pi/90*two_theta/2);
circum = 2*pi*radius;

% Assume pixels approximately lie on circumeference
pixel_angle = pixel_side./circum*360;

%% Compute center of mass
top_dir = 'E:\PureTiRD_full';
dset_name = 'ring%i_zero';
prefix = 'mmpad_img_%i.mat';

R = 4;
N = 546;
COM = zeros(R,N);
angleCOM = zeros(R,N);
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
    angleCOM(ring,:) = atan((radius(ring) + COM(ring,:)-COM(ring,1))/detec_dist);
end


%% Load AWMV
datadir = 'C:\Users\dpqb1\Desktop\mmpad_rd_figures';
fileBase = 'ring%i_zero_couple_fit_data.mat';
outName =  'ring%i_CoM_AWMV.png';
for ring = 1:R
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
datadir = 'C:\Users\dpqb1\Desktop\mmpad_rd_figures';
fileBase = 'ring%i_zero_couple_fit_data.mat';
outName =  'allAWMV.png';
fig_out = figure(5);
for ring = 1:R
    
    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(awmv_az(select_ind,:),'LineWidth',1.5)
    ylabel('AWMV')
    xlabel('time')

end
legend('ring 1','ring 2','ring 3','ring 4','Location','Best')
saveas(fig_out,fullfile(datadir,outName))

%% Plot all CoM
outName =  'allCoM.png';
fig_out = figure(6);
for ring = 1:R
    
    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(COM(ring,:),'LineWidth',1.5)
    ylabel('CoM')
    xlabel('time')

end
legend('ring 1','ring 2','ring 3','ring 4','Location','Best')
saveas(fig_out,fullfile(datadir,outName))

%% Plot all Indep AWMVs
outName =  'indepAWMV.png';
fig_out = figure(7);
for ring = 1:R
    
    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(awmv_az_init,'LineWidth',1.5)
    ylabel('AWMV')
    xlabel('time')

end
legend('ring 1','ring 2','ring 3','ring 4','Location','Best')
saveas(fig_out,fullfile(datadir,outName))
