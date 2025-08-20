function results = compute_metrics(outputDir,sigma,dataset,useMin,indep_only)
%param_select_3D 
if nargin < 5
    indep_only = false;
end
if nargin < 4
    useMin = false;
end

[~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigma,dataset);

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
D_error = zeros(NN,1);
X_error = zeros(NN,1);
vdf_error = zeros(NN,1);
shift_error = zeros(NN,1);

for i = 1:numel(matFileNames)
    % Load outputs
    load(fullfile(outputDir,matFileNames{i}))
    if indep_only
        if outputs.lambda2 > 0
            continue
        end
    end
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
        X = outputs.Y;
    end
    scales = outputs.scales;
    center = (M+1)/2;

    % Compute recons
    AD = reSampleCustomArrayCenter3(N,D,scales,center);
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
        [D_perm, ~] = align_third_dim_and_shift(D, Dtrue);
    
        % Compute errors on recovered X and D 
        D_error(i) = sqrt(sum((D_perm-Dtrue).^2,'all'))/sqrt(sum((Dtrue).^2,'all'));
       
        vdf = sum(X,[1,2]);
        vdf_true = sum(Xtrue,[1,2]);
        vdf_error(i) = sqrt(sum((vdf-vdf_true).^2,'all'))/sqrt(sum((vdf_true).^2,'all'));
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

results = struct();

results.error       = error;
results.rel_error   = rel_error;
results.true_error  = true_error;
results.l1_norm     = l1_norm;
results.l0_norm     = l0_norm;
results.log_penalty = log_penalty;
results.of_penalty  = of_penalty;
results.hs_penalty  = hs_penalty;
results.lambda_vec  = lambda_vec;
results.D_error     = D_error;
results.vdf_error   = vdf_error;
