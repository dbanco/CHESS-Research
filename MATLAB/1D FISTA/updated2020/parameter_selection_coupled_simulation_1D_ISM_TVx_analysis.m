%% Parameter selection
clear all
sM = 1;
eM = 30;
for ijk = [3]
    close all

    % dataset = '/cluster/home/dbanco02/mmpad_polar/ring1_zero/';
    % output_dir = '/cluster/shared/dbanco02/mmpad_1D_indep_param_1/'.;
    
    dset_name = 'gnoise4_nonorm';
    dataset_num = num2str(ijk);
    num_ims = 20;
    home_dir = 'D:\CHESS_data\';
    dset_fit = [dset_name,'_coupled_ISM_TVx5'];
    indep_dir =  'ADMM_Sherman_indep5\';
    dataset =    [home_dir,'simulated_two_spot_1D_',dset_name,'_',dataset_num,'\'];
    init_dir =   [home_dir,indep_dir,'simulated_two_spot_1D_',dset_name,'_',dataset_num,'_simul_init\'];
    output_dir = [home_dir,dset_fit,'\simulated_two_spot_1D_',dset_name,'_',dataset_num,'_coupled_'];

%     datadir = ['/cluster/shared/dbanco02/',dset_name,'_coupled_TV1/'];
%     dataset = ['/cluster/shared/dbanco02/simulated_two_spot_1D_',dset_name,'_',dataset_num,'/'];
%     init_dir =      [datadir,'simulated_two_spot_1D_',dset_name,'_',dataset_num,'_simul_init/'];
%     output_dir =    [datadir,'simulated_two_spot_1D_',dset_name,'_',dataset_num,'_coupled_'];

    baseFileName = 'fista_fit_%i_%i.mat';

    % Universal Parameters
    % Ring sampling parameters
    prefix = 'mmpad_img';
    load([output_dir,'1_final\',sprintf(baseFileName,1,1)])
    polar_image = zeroPad(polar_image,P.params.zeroPad);
    lambda2_vals = Pc.lambda2_values;
    M = numel(lambda2_vals);
    sM = 1;
    eM = M;
    % Construct dictionary
    switch P.basis
        case 'norm2'
            A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
    end
    A0 = unshifted_basis_vector_stack_norm2_zpad(P);

    % Construct distance matrix
    Threshold = 32;
    maxNorm = 0;
    D = constructDistanceMatrix_1D(P,Threshold,maxNorm);
%     maxNorm = 1;
%     D2 = constructDistanceMatrix_1D(P,Threshold,maxNorm);

    % Get error, sparsity, awmv
    err_select = zeros(M,num_ims);
    l0_select = zeros(M,num_ims);
    l1_select = zeros(M,num_ims);
    awmv_az = zeros(M,num_ims);
    x_hats = cell(M,num_ims);
    obj1 = zeros(M,1);
    obj2 = zeros(M,1);
    obj3 = zeros(M,num_ims-1);
    obj3b = zeros(M,num_ims-1);
    
    for k = sM:eM
        fprintf('%i of %i\n',k,M)
        for j = 1:num_ims
            load(fullfile([output_dir,num2str(k),'_final\'],sprintf(baseFileName,1,j)))

            % Fit objective
            fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_hat),P.params.zeroMask);
            polar_vector = polar_image;
            b = polar_vector;
            b = zeroPad(b,P.params.zeroPad);
            x_hats{k,j} = x_hat;

            err_select(k,j) = norm(b(:)-fit(:));
            obj1(k) = obj1(k) + 0.5*norm(b(:)-fit(:))^2;
            
            % Vdf
            az_signal = squeeze(sum(x_hat,1));
            var_sum = sum(az_signal(:));
            vdfs(:,k,j) = az_signal./var_sum;
            awmv_az(k,j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
            
            % Sparsity objective
            l0_select(k,j) = sum(x_hat(:)>0);
            l1_select(k,j) = sum(x_hat(:));
            obj2(k) = obj2(k) + l1_select(k,j)*P.params.lambda;


        end
    end

    % Compute tv objective
    for k = sM:eM
        fprintf('TVx %i of %i\n',k,M)
        for j = 1:num_ims-1
            tv_dist = sum(abs( x_hats{k,j}(:) - x_hats{k,j+1}(:) ));
            obj3(k,j) = tv_dist;
        end
    end

    % Load statistics for independently fit data
    awmv_az_init = zeros(num_ims,1);
    err_indep = zeros(num_ims,1);
    l0_indep = zeros(num_ims,1);
    l1_indep = zeros(num_ims,1);
    vdfs_indep = zeros(P.num_var_t,num_ims);
    tv_indep = zeros(num_ims-1,1);
    x_indep = cell(num_ims,1);
    
    for j = 1:num_ims
        load(fullfile(init_dir,sprintf(baseFileName,1,j)))
        x_indep{j} = x_hat;
        err_indep(j) = err(end-1);
        l0_indep(j) = sum(x_hat(:)>0);
        l1_indep(j) = sum(x_hat(:));
        az_signal = squeeze(sum(x_hat,1));
        var_sum = sum(az_signal(:));
        vdfs_indep(:,j) = az_signal./var_sum;
        awmv_az_init(j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end

    for j = 1:num_ims-1
        tv_dist = sum( abs( x_indep{j}(:) - x_indep{j+1}(:) ));
        tv_indep(j) = tv_dist;
    end

    %% Plot awmv
    figure_dir = ['C:\Users\dpqb1\OneDrive\Desktop\',dset_fit,'_figures\'];
    mkdir(figure_dir)
    % Load truth AWMV
    load([home_dir,'simulated_two_spot_1D_',dset_name,'_11\synth_data.mat'])
    awmv_truth = zeros(num_ims,1);
    for i = 1:num_ims
        awmv_truth(i) = mean(synth_sample{i}.std_theta); 
    end

    close all

    % Plot AWMV
    awmv_fig = figure(1);
    legend_str = {};
    legend_str{1} = 'truth';
    legend_str{2} = '0';
    hold on
    plot(awmv_truth,'LineWidth',1.5)
    plot(awmv_az_init,'LineWidth',1.5)
    kk = 3;
    for k = sM:eM
        hold on
        plot(awmv_az(k,:),'LineWidth',1.5)
        legend_str{kk} = sprintf('%0.03f',lambda2_vals(k));
        kk = kk + 1;
    end
    ylim([0 40])
    ylabel('AWMV_\eta','FontSize',20)
    xlabel('t','FontSize',20)
%     legend(legend_str,'location','best','FontSize',16)
    saveas(awmv_fig,[figure_dir,'awmv_all_',dset_name,'_',dataset_num,'.png'])
    %% Save images

    % Plot err
    mean_err = mean(err_select,2);
    figure(2)
    hold on
    plot(0,mean(err_indep),'o')
    plot(lambda2_vals(sM:eM),mean_err(sM:eM),'o-')
    ylabel('Average Error')
    xlabel('Coupling parameter')

    % Plot l1 norm
    mean_l1 = mean(l1_select,2);
    figure(3)
    hold on
    plot(0,mean(l1_indep),'o')
    plot(lambda2_vals(sM:eM),mean_l1(sM:eM),'o-')
    ylabel('l1-norm')
    xlabel('Coupling parameter')

    % Plot number nonzeros coefficients
    mean_l0 = mean(l0_select,2);
    figure(4)
    loglog(0,mean(l0_indep),'o')
    hold on
    plot(lambda2_vals(sM:eM),mean_l0(sM:eM),'o-')
    ylabel('l0-norm')
    xlabel('Coupling parameter')

    % Plot tv
    tv_time = sum(obj3,1);
    tv_total = sum(obj3,2);
    
    figure(44)
    semilogx(lambda2_vals(sM),sum(tv_indep),'o')
    hold on
    semilogx(lambda2_vals(sM:eM),tv_total(sM:eM),'o-')
    ylabel('TV')
    xlabel('Coupling parameter')

    obj_part1_sum = mean_l1+ 10*mean_err;

    % Plot wass dist over imsages
    figure(441)
    hold on
    plot(tv_indep,'-')
    plot(tv_time,'-')
    legend('indep','regularized')
    ylabel('TV')
    xlabel('time')
%     % Plot AWMV
%     figure(5)
%     legend_str = {};
%     legend_str{1} = 'indep';
%     hold on
%     plot(awmv_az_init,'-','LineWidth',1.5)
%     kk = 2;
%     for k = [4,5,M]
%         hold on
%         plot(awmv_az(sort_i(k),:),'-','LineWidth',1.5)
%         legend_str{kk} = sprintf('%0.04f',sort_lambda2(k));
%         kk = kk + 1;
%     end
%     legend(legend_str,'Location','Best')

    % % Plot L-curve 1
    % figure(6)
    % plot(obj2(sort_i),obj1(sort_i),'o-')
    % xlabel('l1-norm')
    % ylabel('Error')

    % Plot L-curve 2

    Lcurve_fig = figure(7);
    loglog(obj1(sM:eM),sum(obj3(sM:eM,:),2),'o-')
    ylabel('TV')
    xlabel('Error')
    
    %
    miny = min(sum(obj3,2));
    minx = min(obj1);
    coord = [obj1./minx,sum(obj3,2)./miny];
    origin_dist = coord-1;
    [val,select_ind] = min(sum(origin_dist.^2,2));
    
    selectx = obj1(select_ind);
    selecty = sum(obj3(select_ind,:),2);
    hold on
    loglog(selectx,selecty,'s','Markersize',14)
    saveas(Lcurve_fig,[figure_dir,'Lcurve_',dset_name,'_',dataset_num,'.png'])

    % % Plot L-curve 3
    % figure(8)
    % plot(obj2(sort_i),obj3(sort_i),'o-')
    % ylabel('Wasserstein distance')
    % xlabel('l1-norm')
    % 
    % % Plot L-curve 4
    % figure(9)
    % loglog(obj_part1_sum,sum(obj3(sort_i,:),2),'o-')
    % ylabel('Wasserstein distance')
    % xlabel('Error+l1-norm')

    % Plot paramter selected
    select_fig = figure(11);
    legend_str = {};
    legend_str{1} = 'truth';
    legend_str{2} = 'indep';
    hold on
    plot(awmv_truth,'LineWidth',1.5)
    plot(awmv_az_init,'LineWidth',1.5)
    kk = 3;
    for k = [select_ind]
        hold on
        plot(awmv_az(k,:),'LineWidth',1.5)
        legend_str{kk} = sprintf('%0.03f',lambda2_vals(k));
        kk = kk + 1;
    end
    ylim([0 40])
    ylabel('AWMV_\eta','FontSize',20)
    xlabel('t','FontSize',20)
    legend(legend_str,'location','best','FontSize',16)
%     saveas(select_fig,[figure_dir,'awmv_select_',dset_name,'_',dataset_num,'.png'])

    %% Plot vdf surface
    vdf_time_all_fig = figure(56);
    [ha3, pos3] = tight_subplot(2,ceil(M/2)+1,[0.1 0.03],[.02 .08],[.02 .02]); 
    im_ind = 1;
        % Plot surface
    axes(ha3(im_ind))
    imagesc(squeeze(vdfs_indep'))
    shading interp
    caxis([0 0.6])
    colormap(jet)

    title(['\lambda_2 = ','0'])
    ylabel('t')
    xlabel('\sigma')
    
    im_ind = im_ind + 1;
    
    for trial_k = sM:eM
        fprintf('%i of %i\n',trial_k,M)
        for image_num = 1:num_ims

            load(fullfile([output_dir,num2str(trial_k),'a\'],sprintf(baseFileName,1,image_num)))

            polar_vector = squeeze(sum(polar_image,1));
            fit = Ax_ft_1D(A0ft_stack,x_hat)*norm(polar_vector(:));
            az_signal = squeeze(sum(x_hat,1));
            var_sum = sum(az_signal(:));
            vdf_time(trial_k,image_num,:) = az_signal/var_sum;
        end

        % Plot surface
        axes(ha3(im_ind))
        imagesc(squeeze(vdf_time(trial_k,:,:)))
        shading interp
        caxis([0 0.6])
        colormap(jet)

        title(['\lambda_2 = ',sprintf('%1.1d',lambda2_vals(trial_k))])
        ylabel('t')
        xlabel('\sigma')

        im_ind = im_ind + 1;
    end
    saveas(vdf_time_all_fig,[figure_dir,'vdf_time_all_',dset_name,'_',dataset_num,'.png'])

    vdf_time_fig = figure(566)
    trial_k = select_ind;
    fprintf('%i of %i\n',trial_k,M)

    for image_num = 1:num_ims

        load(fullfile([output_dir,num2str(trial_k),'_final\'],sprintf(baseFileName,1,image_num)))

        polar_vector = squeeze(sum(polar_image,1));
        fit = Ax_ft_1D(A0ft_stack,x_hat)*norm(polar_vector(:));
        az_signal = squeeze(sum(x_hat,1));
        var_sum = sum(az_signal(:));
        vdf_time(trial_k,image_num,:) = az_signal/var_sum;
    end
    % Plot surface
    imagesc(squeeze(vdf_time(trial_k,:,:)))
    shading interp
    caxis([0 0.6])
    colormap(jet)

    title(['\lambda_2 = ',sprintf('%1.1d',lambda2_vals(trial_k))])
    ylabel('t')
    xlabel('\sigma')
%     saveas(vdf_time_fig,[figure_dir,'vdf_time_select_',dset_name,'_',dataset_num,'.png'])



    %% Plot fits
    fits_fig = figure(222);
    [ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
    awmv_az_vdfs = zeros(num_ims,1);
    im_ind = 1;
    trial_k = 1;
    for image_num = 1:20

        load(fullfile([output_dir,num2str(trial_k),'_final\'],sprintf(baseFileName,1,image_num)))

        polar_vector = polar_image;
        fit = Ax_ft_1D(A0ft_stack,x_hat);
        az_signal = squeeze(sum(x_hat,1));
        var_sum = sum(az_signal(:));
        awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        b = zeroPad(polar_vector,P.params.zeroPad);

        % Plot
        axes(ha2(im_ind))
        hold on
        plot(b)
        plot(fit)
        legend(sprintf('%i',sum(x_hat(:)>1e-4)),'location','northeast')
        im_ind = im_ind + 1;
    end
%     saveas(fits_fig,[figure_dir,'fits_',dset_name,'_',dataset_num,'.png'])


    %% Plot vdfs
    % figure(333)
    % [ha3, pos3] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
    % awmv_az_vdfs = zeros(num_ims,1);
    % im_ind = 1;
    % trial_k = 5;
    % for image_num = 1:10
    %     
    %     load(fullfile([output_dir,num2str(sort_i(trial_k)),'a\'],sprintf(baseFileName,1,image_num)))
    %     
    %     polar_vector = squeeze(sum(polar_image,1));
    %     fit = Ax_ft_1D(A0ft_stack,x_hat)*norm(polar_vector(:));
    %     az_signal = squeeze(sum(x_hat,1));
    %     var_sum = sum(az_signal(:));
    %     awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    %     % Plot
    %     axes(ha3(im_ind))
    %     plot(az_signal/var_sum)
    %     legend(sprintf('%i',sum(x_hat(:)>1e-6)),'location','northeast')
    %     im_ind = im_ind + 1;
    % end
end

