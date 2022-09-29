function [awmv_indep,awmv_coupled] = getAwmvParamSelect(NN,MM,T,top_dir,sim_name,ps_name)

indepPS_dir = fullfile(top_dir,sim_name,'IndepISM');
coupledPS_dir = fullfile(top_dir,sim_name,'CoupledCGTV');

awmv_indep = zeros(T,NN);
awmv_coupled = zeros(T,NN);

lambda2_inds = zeros(NN,1);

for nn = 1:NN
    % Load param search data
    i_data = load(fullfile(indepPS_dir,[ps_name,'_',num2str(nn)]),...
             'X_indep','P');
    i_indices = i_data.P.indep_select_ind;
    var_theta = i_data.P.var_theta;
    try
        c_data = load(fullfile(coupledPS_dir,[ps_name,'_',num2str(nn)]),...
                 'X_coupled','P');
        c_index = c_data.P.coupled_select_ind;
        X_hat = squeeze(c_data.X_coupled(:,:,c_index,:));
    catch
        c_index = selectParamCoupled(MM,top_dir,sim_name,'CoupledCGTV',ps_name,nn);
        load(fullfile(top_dir,sim_name,'CoupledCGTV',...
            [ps_name,'_',num2str(nn),'_',num2str(c_index),'.mat']),'X_hat');
    end
    lambda2_inds(nn) = c_index;
        
    for t = 1:T
        awmv_indep(t,nn) = computeAWMV_1D(i_data.X_indep(:,:,i_indices(t),t),...
                                          var_theta);
    end
    awmv_coupled(:,nn) = computeAWMV_1D(X_hat,...
                                        var_theta);
end
lambda2_inds
    
