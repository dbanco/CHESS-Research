function objective = eval_trials(topDir,n,sigma,dataset,useMin,num_trials)
%param_select_3D 

[~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigma,dataset);

NN = num_trials;
error = zeros(NN,1);
rel_error = zeros(NN,1);
true_error = zeros(NN,1);
l1_norm = zeros(NN,1);
l0_norm = zeros(NN,1);
log_penalty = zeros(NN,1);
of_penalty = zeros(NN,1);
hs_penalty = zeros(NN,1);
lambda_vec = zeros(NN,3);
D_error = zeros(NN,1);
X_error = zeros(NN,1);
vdf_error = zeros(NN,1);
shift_error = zeros(NN,1);

% Extract the file names and store them in a cell array
for i = 1:num_trials
    inDir = [topDir,'\results_trial_',num2str(i),'_sig_',num2str(n)];
    files = dir(fullfile(inDir, '*.mat'));
    matFileNames = {files.name};
  
    for j = 1:numel(matFileNames)
        % Load outputs
        load(fullfile(inDir,matFileNames{j}))
    
        N = outputs.N;
        M = outputs.M;
        T = outputs.T;
        K = outputs.K;
        KJ = size(Xtrue,3);
        J = KJ/K;
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
    
        % Identify correct ordering and shift of learned dictionary and apply it
        if ~isscalar(Dtrue)
            [D_perm, best_perm, shifts, ~] = align_third_dim_and_shift(D, Dtrue);
            X_perm = apply_perm_to_X(X, J, best_perm, shifts);
        
            % Compute errors on recovered X and D 
            Xerr1 = sqrt(sum((X_perm-Xtrue).^2,'all'))/sqrt(sum((Xtrue).^2,'all'));
            Xerr2 = sqrt(sum((X-Xtrue).^2,'all'))/sqrt(sum((Xtrue).^2,'all'));
            X_error(i) = min(Xerr1,Xerr2);
            D_error(i) = sqrt(sum((D_perm-Dtrue).^2,'all'))/sqrt(sum((Dtrue).^2,'all'));
            if Xerr1 < Xerr2
                vdf = sum(X_perm,[1,2]);
                shift = sum(X_perm,[1,3]);
            else
                vdf = sum(X,[1,2]);
                shift = sum(X,[1,3]);  
            end
            vdf_true = sum(Xtrue,[1,2]);
            vdf_error(i) = sqrt(sum((vdf-vdf_true).^2,'all'))/sqrt(sum((vdf_true).^2,'all'));
            
            shift_true = sum(Xtrue,[1,3]);
            shift_error(i) = sqrt(sum((shift-shift_true).^2,'all'))/sqrt(sum((shift_true).^2,'all'));
        end
    
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
    
    objective = struct();
    objective.error = error;
    objective.rel_error = rel_error;
    objective.true_error = true_error;
    objective.l1_norm = l1_norm;
    objective.l0_norm = l0_norm;
    objective.log_penalty = log_penalty;
    objective.of_penalty = of_penalty;
    objective.hs_penalty = hs_penalty;
    objective.X_error = X_error;
    objective.D_error = D_error;
    objective.vdf_error = vdf_error;
    objective.shift_error = shift_error;
end


end