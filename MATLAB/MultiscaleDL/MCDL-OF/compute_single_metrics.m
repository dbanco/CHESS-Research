function metrics = compute_single_metrics(ADf,y,y_true,D,Dtrue,X,Xtrue,a,K,J,M,opt)

% Compute recons
Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
Yhat = gather(Yhat);

% Compute error metrics
metrics.error = sum((squeeze(y)-Yhat).^2,'all');
metrics.rel_error = metrics.error./sum(squeeze(y).^2,'all');
metrics.true_error = sqrt(sum((y_true-Yhat).^2,'all'));

% Identify correct ordering and shift of learned dictionary and apply it
if ~isscalar(Dtrue)
    [D_perm, ~] = align_third_dim_and_shift(D, Dtrue);

    % Compute errors on recovered X and D 
    metrics.D_error = sqrt(sum((D_perm-Dtrue).^2,'all'))/sqrt(sum((Dtrue).^2,'all'));

    vdf = sum(X,[1,2]);
    vdf_true = sum(Xtrue,[1,2]);
    metrics.vdf_error = sqrt(sum((vdf-vdf_true).^2,'all'))/sqrt(sum((vdf_true).^2,'all'));
end

% Compute log penalty
metrics.log_penalty = sum(vec(log(1 + a.*abs(X))));
% Compute L1-norm
metrics.l1_norm = norm(X(:),1);
% Compute L0-norm
metrics.l0_norm = sum(X(:)>0,'all');
% Compute temporal regularization penalty
metrics.reg_penalty = compute_penalty(X,K,J,opt.regularizer);


% X Error metrics
metrics.x_metric = compute_x_metric(Xtrue,X,K,J);
metrics.x_metric2 = norm(X(:)-Xtrue(:));
metrics.wass_dist = wass_distance(Xtrue,X,K,J);


end