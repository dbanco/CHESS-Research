function results = compute_metrics_reg(outputDir,sigma,dataset,useMin,indep_only)
%param_select_3D 
if nargin < 5
    indep_only = false;
end
if nargin < 4
    useMin = false;
end

[~,~,~,~,~,~,Dtrue] = sim_switch_multiscale_dl(sigma,dataset);

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
reg_penalty = zeros(NN,1);
lambda_vec = zeros(NN,2);
D_error = zeros(NN,1);
x_metric = zeros(NN,1);
x_metric2 = zeros(NN,1);
wass_dist = zeros(NN,1);
vdf_error = zeros(NN,1);

for i = 1:numel(matFileNames)
    % Load outputs
    load(fullfile(outputDir,matFileNames{i}))
    if indep_only
        if outputs.lambda2 > 0
            continue
        end
    end
    
    if useMin
        metrics = outputs.metrics_min;
    else
        metrics = outputs.metrics;
    end

    % Compute error
    error(i) = metrics.error;
    rel_error(i) = metrics.rel_error;
    true_error(i) = metrics.true_error;

    % Identify correct ordering and shift of learned dictionary and apply it
    if ~isscalar(Dtrue)    
        % Compute errors on recovered X and D 
        D_error(i) = metrics.D_error;
        vdf_error(i) = metrics.vdf_error;
    end

    % Compute sparsity penalties
    log_penalty(i) = metrics.log_penalty;
    l1_norm(i) = metrics.l1_norm;
    l0_norm(i) = metrics.l0_norm;

    % Compute regularizer penalties
    reg_penalty(i) = metrics.reg_penalty;

    % Compute coefficient metrics
    x_metric(i) = metrics.x_metric;
    x_metric2(i) = metrics.x_metric2;
    wass_dist(i) = metrics.wass_dist;

    % Get lambda parameters
    lambda_vec(i,1) = outputs.lambda;
    lambda_vec(i,2) = outputs.lambda2;  
end

results = struct();

results.error       = error;
results.rel_error   = rel_error;
results.true_error  = true_error;
results.l1_norm     = l1_norm;
results.l0_norm     = l0_norm;
results.log_penalty = log_penalty;
results.reg_penalty  = reg_penalty;
results.lambda_vec  = lambda_vec;
results.D_error     = D_error;
results.vdf_error   = vdf_error;
results.x_metric = x_metric;
results.x_metric2 = x_metric2;
results.wass_dist = wass_dist;
