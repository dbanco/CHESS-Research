function [awmv_coup,awmv_rse_coup,awmv_indep,awmv_rse_indep] = awmvMCResults(NN,nTrials,T,theta_stds,sim)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

top_dir = 'D:\Simulations';
sim_name = [sim,'PoissonNoise3'];
mc_dir = 'mcTrials';

awmv_coup = zeros(NN,nTrials,T);
awmv_rse_coup = zeros(NN,nTrials);
awmv_indep = zeros(NN,nTrials,T);
awmv_rse_indep = zeros(NN,nTrials);


fprintf('Loading and computing awmv/rse: \n')
lambda2_nn = zeros(NN,1);
for nn = 1:NN
    fprintf('%i, ',nn)
    % Directories for MC trials
    alg_name = 'IndepISM';
    indepPS_dir  = fullfile(top_dir,sim_name,alg_name);
    indepTrials_dir = fullfile(indepPS_dir,[mc_dir,'_',num2str(nn)]);
    alg_name = 'CoupledCGTV';
    coupledPS_dir = fullfile(top_dir,sim_name,alg_name); 
    coupledTrials_dir = fullfile(coupledPS_dir,[mc_dir,'_',num2str(nn)]);
    
    for i = 1:nTrials
        % Load independent and compute awmv rse
        load(fullfile(indepTrials_dir,['trial_',num2str(i)]),'awmv');
        awmv_indep(nn,i,:) = awmv;
        awmv_rse_indep(nn,i) = norm(awmv-theta_stds)/norm(theta_stds);
        
        % Load coupled and compute awmv rse
        load(fullfile(coupledTrials_dir,['trial_',num2str(i)]),'awmv','P');
        awmv_coup(nn,i,:) = awmv;
        awmv_rse_coup(nn,i) = norm(awmv-theta_stds)/norm(theta_stds);
    end
    lambda2_nn(nn) = P.params.lambda2;
end
end

