function [lambda_all,objective] = param_select_3D(outputDir,fig_num,criterion,sigma,y_true,useMin,relax_param)
%param_select_3D 

if nargin < 7
    relax_param = 1.05;
end
if nargin < 6
    useMin = false;
end

% Extract the file names and store them in a cell array
files = dir(fullfile(outputDir, '*.mat'));
matFileNames = {files.name};

NN = numel(matFileNames);
error = zeros(NN,1);
rel_error = zeros(NN,1);
true_error = zeros(NN,1);
l1_norm = zeros(NN,1);
l0_norm = zeros(NN,1);
log_penalty = zeros(NN,1);
of_penalty = zeros(NN,1);
hs_penalty = zeros(NN,1);
lambda_vec = zeros(NN,3);
for i = 1:numel(matFileNames)
    % Load outputs
    load(fullfile(outputDir,matFileNames{i}))

    N = outputs.N;
    M = outputs.M;
    T = outputs.T;
    K = outputs.K;
    y = outputs.y;
    if useMin
        D = outputs.Dmin;
        X = outputs.Ymin;
    else
        D = outputs.D;
        X = outputs.X;
    end
    scales = outputs.scales;
    center = (M+1)/2;

    % Compute recons
    AD = reSampleCustomArrayCenter(N,D,scales,center);
    AD = padarray(AD,[0 M-1 0],0,'post');
    ADf = fft2(AD);
    Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
    Yhat = gather(Yhat);
    % Compute error
    error(i) = sum((squeeze(y)-Yhat).^2,'all');
    rel_error(i) = error(i)./sum(squeeze(y).^2,'all');
    true_error(i) = sqrt(sum((y_true-Yhat).^2,'all'));
    % Compute log penalty
    log_penalty(i) = sum(vec(log(1 + outputs.opt.a.*abs(X))));
    % Compute L1-norm
    l1_norm(i) = norm(X(:),1);
    l0_norm(i) = sum(X(:)>0,'all');
    % Compute OF, HS penalties
    Xpad = padarray(squeeze(outputs.X),[1 1 1],0,'pre');
    Fx = diffxHS(Xpad);
    Fy = diffyHS(Xpad);
    Ft = difftHS(Xpad);
    [Jof, Jhs] = HSobjectivePaper(Fx,Fy,Ft,outputs.Uvel,outputs.Vvel,K,outputs.opt.Smoothness/outputs.lambda2);
    of_penalty(i) = Jof;
    hs_penalty(i) = Jhs;
    % Get lambda parameters
    lambda_vec(i,1) = outputs.lambda;
    lambda_vec(i,2) = outputs.lambda2;
    lambda_vec(i,3) = outputs.opt.Smoothness;
end

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
        crit1 = error/(N*T) < relax_param*sigma^2;
        crit2 = error/(N*T) > (2-relax_param)*sigma^2;
        include = crit1 & crit2;
        
        if sum(include) == 0 % default to relaxed discrepancy
            crit = abs(error/(N*T) - relax_param*sigma^2);
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
if fig_num > 0
    figure(fig_num)
    loglog(l1_norm(1:end),error(1:end),'o')
    ylabel('Error')
    xlabel('l_1-norm')
    hold on
    loglog(l1_norm(selInd),error(selInd),'sr','MarkerSize',10)
end

objective = struct();
objective.error = error(i);
objective.rel_error = rel_error(i);
objective.true_error = true_error(i);
objective.l1_norm = l1_norm(i);
objective.log_penalty = log_penalty(i);
objective.of_penalty = of_penalty(i);
objective.hs_penalty = hs_penalty(i);


end