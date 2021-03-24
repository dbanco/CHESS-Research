% MMPAD metadata
mmpad_dims = [396,265];
rings = {'020','112','021','004'};
d_space = [1.27773,1.24845,1.23272,1.17138];
two_theta = [6.89132,7.05316,7.14328,7.5179];

total_time = 136.25;
fps = 4;
x_time = 0:1/fps:total_time;

pixel_side = 150e-6;
gap_width = 0.75e-3;

detec_dist = 4.6;
detect_angle = 14.6;

% Compute circumeference of rings
radius = 4.6*tan(pi/90*two_theta/2);
circum = 2*pi*radius;

% Assume pixels approximately lie on circumeference
pixel_angle = pixel_side./circum*360;


top_dirs = {'D:\MMPAD_data_nr1','D:\MMPAD_data_nr2'};
data_tags = {'_nr1','_nr2'};
num_rings = [4,3];
for d_num = 1
    top_dir = top_dirs{d_num};
    data_tag = data_tags{d_num};
    dset_name = 'ring%i_zero';
    prefix = 'mmpad_img_%i.mat';
    
    %% Compute center of mass
    R = num_rings(d_num);
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
    close all

    %% Plot all CoM and convert to angle units
    outName =  ['allCoM',data_tag,'.png'];
    datadir = 'C:\Users\dpqb1\Desktop\paper_figures';
    fileBase = ['ring%i_zero_couple_fit_data',data_tag,'.mat'];
    fig_CoM = figure(9);
    
    for ring = 1:R
        load(fullfile(datadir,sprintf(fileBase,ring)))
        hold on
        plot(x_time,angleCOM(ring,:),'LineWidth',1.5)
        ylabel('\Delta CoM in 2\theta (degrees)')
        xlabel('time (s)')

    end
    c_leg = legend(rings,'Location','Best');
    ctitle = get(c_leg,'Title');
    set(ctitle,'String','Ring indices \{hkl\}')
    grid on
    saveas(fig_CoM,fullfile(datadir,outName))

    %% Load AWMV
    fileBase = ['ring%i_zero_couple_fit_data',data_tag,'.mat'];
    outName =  ['ring%i_CoM_AWMV', data_tag, '.png'];
    for ring = 1:R
        fig_out = figure(ring);
        load(fullfile(datadir,sprintf(fileBase,ring)))

        % Rescale center of Mass
        minAWMV = min(awmv_az(select_ind,:)*pixel_angle(ring));
        maxAWMV = max(awmv_az(select_ind,:)*pixel_angle(ring));
        scaleCOM = (angleCOM(ring,:)-min(angleCOM(ring,:)))/(max(angleCOM(ring,:))-min(angleCOM(ring,:)));
        scaleCOM = (scaleCOM + minAWMV/maxAWMV)*maxAWMV;

        hold on
        xlabel('time (s)')
        ylabel('AWMV in \eta (degrees)')
        title(sprintf('Ring %i',ring))
        plot(x_time,awmv_az(select_ind,:)*pixel_angle(ring),'LineWidth',1.5)
        plot(x_time,scaleCOM,'LineWidth',1.5)
        legend('AWMV','scaled \Delta CoM','Location','Best')
        grid on
        saveas(fig_out,fullfile(datadir,sprintf(outName,ring)))
    end

    %% PLot all awmv
    fileBase = ['ring%i_zero_couple_fit_data',data_tag,'.mat'];
    outName =  ['allAWMV',data_tag,'.png'];
    fig_out = figure(5);
    for ring = 1:R

        load(fullfile(datadir,sprintf(fileBase,ring)))
        hold on
        plot(x_time,awmv_az(select_ind,:)*pixel_angle(ring),'LineWidth',1.5)
        ylabel('AWMV in \eta (degrees)')
        xlabel('time (s)')

    end
    c_leg = legend(rings,'Location','Best');
    ctitle = get(c_leg,'Title');
    set(ctitle,'String','Ring indices \{hkl\}')
    grid on
    saveas(fig_out,fullfile(datadir,outName))


    %% Plot all Indep AWMVs
    outName = ['indepAWMV',data_tag,'.png'];
    fig_out = figure(7);
    for ring = 1:R

        load(fullfile(datadir,sprintf(fileBase,ring)))
        hold on
        plot(x_time,awmv_az_init*pixel_angle(ring),'LineWidth',1.5)
        ylabel('AWMV in \eta (degrees)')
        xlabel('time (s)')

    end
%     c_leg = legend(rings,'Location','Best');

    grid on
    saveas(fig_out,fullfile(datadir,outName))
    
    figure(7)   
    load('masked_ring4_awmv')
    plot(x_time,masked_ring4_awmv*pixel_angle(4),'LineWidth',1.5)
    new_rings = {rings{1},rings{2},rings{3},rings{4},'004 updated'};
    c_leg = legend(new_rings,'Location','Best');
    ctitle = get(c_leg,'Title');
    set(ctitle,'String','Ring indices \{hkl\}')

%     close all
end
