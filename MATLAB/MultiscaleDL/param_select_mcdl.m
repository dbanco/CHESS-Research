function [lambda_all,objective] = param_select_mcdl(results,criterion,sigma,dataset,relax_param,fig_num,indep_only)
%param_select_3D 
if nargin < 5
    relax_param = 1.05;
end

[~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigma,dataset);

error       = results.error;
rel_error   = results.rel_error;
true_error  = results.true_error;
l1_norm     = results.l1_norm;
l0_norm     = results.l0_norm;
log_penalty = results.log_penalty;
of_penalty  = results.of_penalty;
hs_penalty  = results.hs_penalty;
lambda_vec  = results.lambda_vec;
D_error     = results.D_error;
X_error     = results.X_error;
vdf_error   = results.vdf_error;
shift_error = results.shift_error;

% Normalized origin distance criterion
switch criterion
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
    case 'discrepancy range2'
        if sigma == 0
            crit1 = error/(N*T) < relax_param^2;
            crit2 = crit1;
        else
            crit1 = error/(N*T) < relax_param^2+sigma^2;
            crit2 = error/(N*T) > -relax_param^2+sigma^2;
        end
        include = crit1 & crit2;

        if sum(include) == 0 % default to relaxed discrepancy
            crit = abs(error/(N*T) - relax_param^2+sigma^2);
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
end

if isempty(lambda_all)
    error('No parameter selected')
end

if fig_num > 0
    figure(fig_num)
    loglog(log_penalty(1:end),error(1:end),'o')
    ylabel('Error')
    xlabel('log penalty')
    hold on
    loglog(log_penalty(selInd),error(selInd),'sr','MarkerSize',10)
end

objective = struct();
objective.error = error(selInd);
objective.rel_error = rel_error(selInd);
objective.true_error = true_error(selInd);
objective.l1_norm = l1_norm(selInd);
objective.l0_norm = l0_norm(selInd);
objective.log_penalty = log_penalty(selInd);
objective.of_penalty = of_penalty(selInd);
objective.hs_penalty = hs_penalty(selInd);
objective.X_error = X_error(selInd);
objective.D_error = D_error(selInd);
objective.vdf_error = vdf_error(selInd);
objective.shift_error = shift_error(selInd);


end