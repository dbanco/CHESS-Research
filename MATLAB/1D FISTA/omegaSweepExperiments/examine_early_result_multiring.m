close all
% clear all

% plot indep awmv
% load('E:\MMPAD_omega\coupleld4_ex_data\awmv_indep_ring1omega5')
% plot(awmv_indep(:,40))

resultName = '_coupled_CG_TVphi_Mirror11\coupled_fit_28.mat';

[pixel_angle, x_time] = angleConversionMMPAD();

delta_awmv_all = zeros(546,4,4);
for om = 2:5
    for ring = 1:4
        if om == 5
            om5_ind = [29 23 25 25];
            load(fullfile('C:\Users\dpqb1\Desktop\AWMV_mirror_Figures',...
                 ['ring',num2str(ring),'_zero_mirror_coupled_awmv']),'awmv_az')
            awmv = awmv_az(om5_ind(ring),:);
        else
            topDir = 'E:\MMPAD_omega'; %omD = 
            load(['E:\MMPAD_omega\coupleld4_ex_data\ring',...
                   num2str(ring),'omega',num2str(om),resultName])
            dirSplit = strsplit(P.dataset,'/');
            dataDir = fullfile(topDir,dirSplit{end-1},dirSplit{end});
            awmv = computeAWMV_1D(X_hat,P.var_theta);
        end
        delta_awmv_all(:,ring,om-1) = (awmv - min(awmv))*pixel_angle(ring);
        
        % To looke at original coupled fit instead
        % topDir = 'E:\MMPAD_data_nr1';%ogD = 
        % load('E:\MMPAD_data_nr1\ring1_zero_coupled_CG_TVphi_Mirror5\coupled_fit_29.mat')
        % dirSplit = strsplit(P.dataset,'/');
        % dataDir = fullfile(topDir,dirSplit{end});
        % flipit = 1;

        % Plot Fits
%         T = P.num_ims;
%         Time = 1:15:120;
%         B = loadMMPADmirror(dataDir,Time,P);
%         A0 = unshifted_basis_vector_ft_stack_zpad(P);
%         fig3 = plotSomeFits(X_hat,B,A0,Time);


    end
end

%% Plot AWMV by Omega
close all
savePNG = 0;
plotDeltaAWMVbyOmega(x_time,delta_awmv_all,savePNG);
%% Plot AWMV by Ring
close all
savePNG = 0;
plotDeltaAWMVbyRing(x_time,delta_awmv_all,savePNG);
%% Plot Error
figure
i = 1;
T = 546;
A0 = unshifted_basis_vector_ft_stack_zpad(P);
for om = 2:5
    for ring = 1:4
       
        if om == 5
            om5_ind = [29 23 25 25];
            load(fullfile(['E:\MMPAD_data_nr1\ring',num2str(ring),...
                 '_zero_coupled_CG_TVphi_Mirror5'],...
                 ['coupled_fit_',num2str(om5_ind(ring))]))
        else
            load(['E:\MMPAD_omega\coupleld4_ex_data\ring',...
                num2str(ring),'omega',num2str(om),resultName])
        end
        dataDir = fullfile(topDir,['omega',num2str(om)],...
                                  ['ring',num2str(ring)]);
        B = loadMMPADmirror(dataDir,1:T,P);
        if om == 5
            B = flip(B);
        end
        relErr = computeRelErr(B,X_hat,A0);
        
        subplot(4,4,i)
        plot(relErr)
        title('Relative Error')
        i = i + 1;
    end
end


