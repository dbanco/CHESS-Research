function [objectives_of,objectives_indep] = compute_objectives_trials(sig_ind,sigmas,dataset,useMin,num_trials,topDir,topDir2,selected_lam_s,selected_lam_of,selected_lam_hs,lambdaVals,lambdaOFVals,lambdaHSVals)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

objectives_indep = cell(numel(sig_ind),1);
objectives_of = cell(numel(sig_ind),1);
for n = sig_ind

    j_s_select = selected_lam_s(n);
    j_of_select = selected_lam_of(n);
    j_hs_select = selected_lam_hs(n);
    
    lambda_all = [lambdaVals(j_s_select),...
                  lambdaOFVals(j_of_select),...
                  lambdaHSVals(j_hs_select)];
    HSiters = 100;

    fprintf('sig_ind = %i, ',n)

    fprintf('Indep, ')
    objFile = [topDir,'\obj_indep_sig_',num2str(n),'.mat'];
    if exist(objFile,'file')
        load(objFile)
    else    
        lambda_inds = [j_s_select,1,1];
        objective_indep = eval_trials(topDir,n,sigmas(n),dataset,useMin,num_trials,...
            true,lambda_inds,lambda_all,HSiters);
        objectives_indep{n} = objective_indep;  
        save(objFile,'objectives_indep')
    end
    
    fprintf('OF\n')
    objFile = [topDir2,'\obj_of_sig_',num2str(n),'.mat'];
    if exist(objFile,'file')
        load(objFile)
    else    
        lambda_inds = [j_s_select,j_of_select,j_hs_select];
        objective_of = eval_trials(topDir2,n,sigmas(n),dataset,useMin,num_trials,...
            false,lambda_inds,lambda_all,HSiters);
        objectives_of{n} = objective_of;  
        save(objFile,'objectives_of')
    end
end

end