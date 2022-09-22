function [awmv_indep,awmv_coupled] = getAwmvParamSelect(sim,NN,T)

top_dir = 'D:\Simulations';
sim_name = [sim,'PoissonNoise3'];
ps_name = 'lineSearchNoise';
indepPS_dir = fullfile(top_dir,sim_name,'IndepISM');
coupledPS_dir = fullfile(top_dir,sim_name,'CoupledCGTV');

awmv_indep = zeros(T,NN);
awmv_coupled = zeros(T,NN);

for nn = 1:NN
    % Load param search data
    i_data = load(fullfile(indepPS_dir,[ps_name,'_',num2str(nn)]),...
             'X_indep','P');
    c_data = load(fullfile(coupledPS_dir,[ps_name,'_',num2str(nn)]),...
             'X_coupled','P');
    i_indices = i_data.P.indep_select_ind;
    c_index = c_data.P.coupled_select_ind;
    for t = 1:T
        awmv_indep(t,nn) = computeAWMV_1D(i_data.X_indep(:,:,i_indices(t),t),...
                                 i_data.P.var_theta);
    end
    awmv_coupled(:,nn) = computeAWMV_1D(squeeze(c_data.X_coupled(:,:,c_index,:)),...
                                             c_data.P.var_theta);
end
    
