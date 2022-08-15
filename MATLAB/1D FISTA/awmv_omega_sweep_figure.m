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

R = 4;

%% AWMV indep
datadir = 'E:\MMPAD_omega';
fileBase = 'ring%iomega%i_mirror_indep_awmv.mat';
outName =  'ring%iomega%i_indep_AWMV.png';
indepBase = 'ring%i_zero_mirror_indep_awmv.mat';
for ring = 1:R
    fig_out = figure(ring);
    hold on
    for om = 2:4
        load(fullfile(datadir,sprintf(fileBase,ring,om)))
        % Rescale center of Mass
        minAWMV = min(awmv);
        maxAWMV = max(awmv);
        plot(awmv,'LineWidth',1.5)
    end
    load(fullfile(datadir,sprintf(indepBase,ring)))
    plot(awmv,'LineWidth',1.5)
    title(sprintf('Ring %i',ring))
    legend('2^\circ','3^\circ','4^\circ','5^\circ')
%     saveas(fig_out,fullfile(datadir,sprintf(outName,ring,om)))
end

%% AWMV coupled
datadir = 'E:\MMPAD_omega\coupled';
fileBase = 'ring%iomega%i_mirror_coupled_awmv.mat';
outName =  'ring%iomega%i_coupled_AWMV.png';
origdir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
origBase = 'ring%i_zero_mirror_coupled_awmv.mat';
orig_ind = [29 23 25 25];
for ring = 1:R
    fig_out = figure(ring);
    hold on
    for om = 2:4
        load(fullfile(datadir,sprintf(fileBase,ring,om)))
        % Rescale center of Mass
        awmv = awmv_az(select_ind,:);
        minAWMV = min(awmv);
        maxAWMV = max(awmv);
        plot(awmv,'LineWidth',1.5)
    end
    load(fullfile(origdir,sprintf(origBase,ring)))
    plot(awmv_az(orig_ind(ring),:),'LineWidth',1.5)
    title(sprintf('Ring %i',ring))
    legend('2^\circ','3^\circ','4^\circ','5^\circ')
    saveas(fig_out,fullfile(datadir,sprintf(outName,ring,om)))
end

%% Compare fits
datadir = 'E:\MMPAD_omega\coupled';
fileBase = 'ring%iomega%i_mirror_coupled_awmv.mat';
outName =  'ring%iomega%i_coupled_AWMV.png';
origdir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
origBase = 'ring%i_zero_mirror_coupled_awmv.mat';

om_dir = {'tt','omega2','omega3','omega4'};
r_dir = {'ring1','ring2','ring3','ring4'};

orig_ind = [29 23 25 25];
for ring = 2
    fig_out = figure(ring);
    hold on
    for om = 2:4
        load(fullfile(datadir,sprintf(fileBase,ring,om)))   
        
        figure(ring)
        % Rescale center of Mass
        awmv = awmv_az(select_ind,:);
        minAWMV = min(awmv);
        maxAWMV = max(awmv);
        plot(awmv,'LineWidth',1.5)
        
        % Plot fits
        output_name = '_coupled_CG_TVphi_Mirror';
        output_subdir = [r_dir{ring},om_dir{om},output_name];
        output_dir  = fullfile(datadir,output_subdir);
        baseFileName = 'coupled_fit_%i.mat';
        load(fullfile(output_dir,sprintf(baseFileName,1)))
        [N,K,T] = size(X_hat);
        B = zeros(N,T);
        dataset = fullfile('E:\MMPAD_omega',om_dir{om},r_dir{ring});
        
        tt = 28;
        for j = tt
            b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
            b = P.dataScale*sum(b_data.polar_image,2);

            % Mirror data
            nn = numel(b);
            pad1 = floor(nn/2);
            pad2 = ceil(nn/2);
            N = nn + pad1 + pad2;
            b_mirror = zeros(N,1);
            b_mirror((pad1+1):(pad1+nn)) = b;
            b_mirror((1+N-pad2):N) = flip(b((nn-pad2+1):nn));
            b_mirror(1:pad1) = flip(b(1:pad1));
            B(:,j) = b_mirror;
        end
        % Load coupled fit for selected ind
        e_data = load(fullfile(output_dir,sprintf(baseFileName,select_ind)),'P','X_hat');
        A0ft_stack = unshifted_basis_vector_ft_stack_zpad(e_data.P);
        tt = 200;
        b = B(:,tt);
        x = squeeze(e_data.X_hat(:,:,tt));
        fit = Ax_ft_1D(A0ft_stack,x);
        figure(om+9)
        hold on
        plot(b)
        plot(fit)
    end
    load(fullfile(origdir,sprintf(origBase,ring)))
    figure(ring)
    plot(awmv_az(orig_ind(ring),:),'LineWidth',1.5)
    title(sprintf('Ring %i',ring))
    legend('2^\circ','3^\circ','4^\circ','5^\circ')
%     saveas(fig_out,fullfile(datadir,sprintf(outName,ring,om)))
end

%% AWMV coupled
datadir = 'E:\MMPAD_omega\coupled';
fileBase = 'ring%iomega%i_mirror_coupled_awmv.mat';
outName =  'ring%iomega%i_coupled_AWMV.png';
origdir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
origBase = 'ring%i_zero_mirror_coupled_awmv.mat';
orig_ind = [29 23 25 25];
for ring = 2
    fig_out = figure(ring);
    hold on
    for om = 4
        load(fullfile(datadir,sprintf(fileBase,ring,om)))
        % Rescale center of Mass
        awmv = awmv_az(30,:);
        minAWMV = min(awmv);
        maxAWMV = max(awmv);
        plot(awmv,'LineWidth',1.5)
    end
    load(fullfile(origdir,sprintf(origBase,ring)))
    plot(awmv_az(orig_ind(ring),:),'LineWidth',1.5)
    title(sprintf('Ring %i',ring))
    legend('2^\circ','3^\circ','4^\circ','5^\circ')
%     saveas(fig_out,fullfile(datadir,sprintf(outName,ring,om)))
end


%% Plot fit from orig and smaller sweep
datadir = 'E:\MMPAD_omega\coupled';
fileBase = 'ring%iomega%i_mirror_coupled_awmv.mat';
outName =  'ring%iomega%i_coupled_AWMV.png';
origdir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
indepBase = 'ring%i_zero_mirror_coupled_awmv.mat';

om_dir = {'tt','omega2','omega3','omega4'};
r_dir = {'ring1','ring2','ring3','ring4'};

orig_ind = [29 23 25 25];
for ring = 2
    fig_out = figure(ring);
    hold on
    for om = 2
        load(fullfile(datadir,sprintf(fileBase,ring,om)))   
        
        figure(ring)
        % Rescale center of Mass
        awmv = awmv_az(select_ind,:);
        minAWMV = min(awmv);
        maxAWMV = max(awmv);
        plot(awmv,'LineWidth',1.5)
        
        % Plot fits
        output_name = '_coupled_CG_TVphi_Mirror';
        output_subdir = [r_dir{ring},om_dir{om},output_name];
        output_dir  = fullfile(datadir,output_subdir);
        baseFileName = 'coupled_fit_%i.mat';
        load(fullfile(output_dir,sprintf(baseFileName,1)))
        [N,K,T] = size(X_hat);
        B = zeros(N,T);
        dataset = fullfile('E:\MMPAD_omega',om_dir{om},r_dir{ring});
        
        tt = 28;
        for j = tt
            b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
            b = P.dataScale*sum(b_data.polar_image,2);
            B(:,j) = mirrorData(b);
        end
        % Load coupled fit for selected ind
        e_data = load(fullfile(output_dir,sprintf(baseFileName,select_ind)),'P','X_hat');
        A0ft_stack = unshifted_basis_vector_ft_stack_zpad(e_data.P);
        tt = 200;
        b = B(:,tt);
        x = squeeze(e_data.X_hat(:,:,tt));
        fit = Ax_ft_1D(A0ft_stack,x);
        figure(om+9)
        hold on
        plot(b)
        plot(fit)
    end
    load(fullfile(origdir,sprintf(indepBase,ring)))
    figure(ring)
    plot(awmv_az(orig_ind(ring),:),'LineWidth',1.5)
    title(sprintf('Ring %i',ring))
    legend('2^\circ','3^\circ','4^\circ','5^\circ')
%     saveas(fig_out,fullfile(datadir,sprintf(outName,ring,om)))
end

%% Plot AWMV for different parameter values
datadir = 'E:\MMPAD_omega\coupled';
fileBase = 'ring%iomega%i_mirror_coupled_awmv.mat';
outName =  'ring%iomega%i_coupled_AWMV.png';
origdir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
indepBase = 'ring%i_zero_mirror_coupled_awmv.mat';
orig_ind = [29 23 25 25];
for ring = 1
    fig_out = figure(ring);
    hold on
    for om = 4
        load(fullfile(datadir,sprintf(fileBase,ring,om)))
        % Rescale center of Mass
        awmv = awmv_az(30,:);
        minAWMV = min(awmv);
        maxAWMV = max(awmv);
        plot(awmv,'LineWidth',1.5)
    end
%     load(fullfile(origdir,sprintf(origBase,ring)))
%     plot(awmv_az(orig_ind(ring),:),'LineWidth',1.5)
%     title(sprintf('Ring %i',ring))
%     legend('2^\circ','3^\circ','4^\circ','5^\circ')
%     saveas(fig_out,fullfile(datadir,sprintf(outName,ring,om)))
end

%% Plot fits for different parameter values

datadir = 'E:\MMPAD_omega\coupled';
fileBase = 'ring%iomega%i_mirror_coupled_awmv.mat';
outName =  'ring%iomega%i_coupled_AWMV.png';
origdir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
indepBase = 'ring%i_zero_mirror_coupled_awmv.mat';
orig_ind = [29 23 25 25];

ring = 2;
close all
dataset =  fullfile('E:\MMPAD_data_nr1',sprintf('ring%i_zero',ring));
m = 10;
tt = 200;


omegDir = 'ring%iomega%i_coupled_CG_TVphi_Mirror';
fName = 'coupled_fit_%i.mat';

% Load original fit
fitOrig = fullfile('E:\MMPAD_data_nr1',...
                   'ring2_zero_coupled_CG_TVphi_Mirror7');
ogData = load(fullfile(fitOrig,sprintf(fName,m)));
A0 = unshifted_basis_vector_ft_stack_zpad(ogData.P);

fitog = Ax_ft_1D(A0,ogData.X_hat(:,:,tt));
for j = tt
    b_data = load(fullfile(dataset,[ogData.P.prefix,'_',num2str(j),'.mat']));
    b = ogData.P.dataScale*sum(b_data.polar_image,1);
    b_m = mirrorData(b);
end

figure(1)
subplot(5,1,1)
hold on
plot(b_m)
plot(fitog)
legend('data','og')

subplot(5,1,5)
hold on
plot(b_m)

for om = 2:4
    % Load an omega fit
    fitOmeg = fullfile(datadir,sprintf(omegDir,ring,om));
    omData = load(fullfile(fitOmeg,sprintf(fName,m)));
    dataset =  fullfile('E:\MMPAD_omega',sprintf('omega%i',om),...
                sprintf('ring%i',ring));
    bom_data = load(fullfile(dataset,[ogData.P.prefix,'_',num2str(tt),'.mat']));
    bom = omData.P.dataScale*sum(bom_data.polar_image,2);
    bom_m = mirrorData(bom);
    fitom = Ax_ft_1D(A0,omData.X_hat(:,:,tt));
    subplot(5,1,om)
    hold on
    plot(bom_m)
    plot(fitom)
    legend('data',num2str(om))
    subplot(5,1,5)
    plot(bom_m)
end




%     load(fullfile(origdir,sprintf(origBase,ring)))
%     plot(awmv_az(orig_ind(ring),:),'LineWidth',1.5)
%     title(sprintf('Ring %i',ring))
%     legend('2^\circ','3^\circ','4^\circ','5^\circ')
%     saveas(fig_out,fullfile(datadir,sprintf(outName,ring,om)))


%{
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
%}