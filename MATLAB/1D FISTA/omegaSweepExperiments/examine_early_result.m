% close all
% clear all

topDir = 'E:\MMPAD_omega'; %omD = 
load('E:\MMPAD_omega\coupleld4_ex_data\ring1omega4_coupled_CG_TVphi_Mirror6\coupled_fit_27.mat')
dirSplit = strsplit(P.dataset,'/');
dataDir = fullfile(topDir,dirSplit{end-1},dirSplit{end});
flipit = 0;

% topDir = 'E:\MMPAD_data_nr1';%ogD = 
% load('E:\MMPAD_data_nr1\ring1_zero_coupled_CG_TVphi_Mirror5\coupled_fit_29.mat')
% dirSplit = strsplit(P.dataset,'/');
% dataDir = fullfile(topDir,dirSplit{end});
% flipit = 1;

%% Load corresponding MMPAD data

T = P.num_ims;
Time = [65:128];
TT = numel(Time);
B = loadMMPADmirror(dataDir,Time,P);
A0 = unshifted_basis_vector_ft_stack_zpad(P);

% Plot some fits
figure
i = 1;
for t = Time
    
    b = B(:,t);
    if flipit
        b = flip(b);
    end
    x = squeeze(X_hat(:,:,t));
    fit = Ax_ft_1D(A0,x); 
    
    subplot(ceil(sqrt(TT)),ceil(sqrt(TT)),i)
    hold on
    plot(b)
    plot(fit)
    i = i + 1;
    title(num2str(norm(b-fit)/norm(b)))
end

%% Plot AWMV
figure
hold on
% plot indep awmv
load('E:\MMPAD_omega\coupleld4_ex_data\awmv_indep_ring1omega5')
plot(awmv_indep(:,40))
awmv = computeAWMV_1D(X_hat,P.var_theta);
plot(awmv)
title([dirSplit{end},dirSplit{end-1}])


%% Plot Error
figure
B = loadMMPADmirror(dataDir,1:T,P);
if flipit
    B = flip(B,1);
end
relErr = computeRelErr(B,X_hat,A0);
plot(relErr)
title('Relative Error')

