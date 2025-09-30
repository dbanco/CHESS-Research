function [lambda_all,objective] = param_select_mcdl_reg(results,criterion,sigma,dataset,relax_param,fig_num)
%param_select_3D 
if nargin < 5
    relax_param = 1.05;
end

[~,~,N,M,T,~,~] = sim_switch_multiscale_dl(sigma,dataset);

error       = results.error;
rel_error   = results.rel_error;
true_error  = results.true_error;
l1_norm     = results.l1_norm;
l0_norm     = results.l0_norm;
log_penalty = results.log_penalty;
reg_penalty  = results.reg_penalty;
lambda_vec  = results.lambda_vec;
D_error     = results.D_error;
vdf_error   = results.vdf_error;
x_metric = results.x_metric;
x_metric2 = results.x_metric2;
wass_dist = results.wass_dist;

% Normalized origin distance criterion
switch criterion
    case 'triangle'
        selInd = find_Lcurve_kink_triangle(error,log_penalty,lambda_vec(:,1));
        lambda_all = lambda_vec(selInd,:);
    case 'l-curve'
        selInd = find_Lcurve_kink(error,log_penalty,lambda_vec(:,1));
        lambda_all = lambda_vec(selInd,:);
    case 'discrepancy'
        crit = abs(error/(N*T) - sigma^2);
        [~,selInd] = min(crit);
        lambda_all = lambda_vec(selInd,:);
    case 'relaxed discrepancy'
        crit = abs(error/(N*T) - relax_param*sigma^2);
        [~,selInd] = min(crit);
        lambda_all = lambda_vec(selInd,:);
    case 'discrepancy range'
        if sigma == 0
            crit1 = error/(N*T) < 0.0025^2;
            crit2 = crit1;
        else
            crit1 = error/(N*T) < relax_param*sigma^2;
            crit2 = error/(N*T) > (2-relax_param)*sigma^2;
        end
        include = crit1 & crit2;

        if sum(include) == 0 % default to relaxed discrepancy
            crit = abs(error/(N*T) - relax_param*sigma^2);
            [~,selInd] = min(crit);
            lambda_all = lambda_vec(selInd,:);
            % error('No solution in discrepancy range');
        else 
            % crit3 = l0_norm;
            crit3 = log_penalty;
            exclude = ~include;
            crit3(exclude) = (N+M-1)*T;
            [~,selInd] = min(crit3);
            lambda_all = lambda_vec(selInd,:);
        end
    case 'discrepancy range mmpad'
        mmpad_thresh = 1/(N*T);
        crit1 = error/(N*T) < relax_param*mmpad_thresh;
        crit2 = error/(N*T) > (2-relax_param)*mmpad_thresh;
        include = crit1 & crit2;
        
        if sum(include) == 0 % default to relaxed discrepancy
            crit = abs(error/(N*T) - relax_param*mmpad_thresh);
            [~,selInd] = min(crit);
            lambda_all = lambda_vec(selInd,:);
            error('No solution in discrepancy range');
        else 
            crit3 = l0_norm;
            exclude = ~include;
            crit3(exclude) = numel(X);
            [~,selInd] = min(crit3);
            lambda_all = lambda_vec(selInd,:);
        end

    case 'truth_error'
        [~,selInd] = min(true_error);
        lambda_all = lambda_vec(selInd,:);
    case 'x_metric'
        [~,selInd] = min(x_metric);
        lambda_all = lambda_vec(selInd,:);
    case 'x_metric2'
        [~,selInd] = min(x_metric2);
        lambda_all = lambda_vec(selInd,:);
    case 'wass_dist'
        [~,selInd] = min(wass_dist);
        lambda_all = lambda_vec(selInd,:);    
end

if isempty(lambda_all)
    error('No parameter selected')
end

if fig_num > 0
    figure(fig_num)
    loglog(reg_penalty(1:end),error(1:end),'o')
    ylabel('Error')
    xlabel('of penalty')
    hold on
    loglog(reg_penalty(selInd),error(selInd),'sr','MarkerSize',10)
end

objective = struct();
objective.error = error(selInd);
objective.rel_error = rel_error(selInd);
objective.true_error = true_error(selInd);
objective.l1_norm = l1_norm(selInd);
objective.l0_norm = l0_norm(selInd);
objective.log_penalty = log_penalty(selInd);
objective.reg_penalty = reg_penalty(selInd);
objective.D_error = D_error(selInd);
objective.vdf_error = vdf_error(selInd);

end