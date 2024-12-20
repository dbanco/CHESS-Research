lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
sigmas = 0:0.01:0.1;
NN = numel(sigmas);

tradeoff_s = 0.4;
tradeoff_of = 1;
scaleP = [0.4,5.4,9.48,116,0,100];

criterion = 'discrepancy';
% criterion = 'truth_error';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);

% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_10_Dtrue_Xtrue';
% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_31_Dtrue1_Xzeros0';
% dirStartS = 'steps_matched_results';
fig_num = 22;

dataset = 'steps_matched';
penalty = 'log';
testType = 'Dflat0_Xzeros0';
topDir = 'E:\Outputs_sim1_trials\';

obj_array = zeros(NN,10);
true_error = zeros(NN,10);
data_error = zeros(NN,10);
l1_norm = zeros(NN,10);
l0_norm = zeros(NN,10);
log_penalty = zeros(NN,10);

for r = 1:10
    for nn = 2:9
        [~,y_true,~,~,~] = gaus_example_switch_multiscale_dl(sigmas(nn),dataset);
        % File
        dir1 = sprintf(['Outputs%i_',testType],r);
        dir2 = sprintf([dataset,'_',penalty,'_results_sig_%i'],nn);
        dir12 = fullfile(topDir,dir1,dir2);
        files = dir(fullfile(dir12, '*.mat'));
        fname = files.name;
        % Load "outputs"
        load(fullfile(dir12,fname));
        [Jfn,Jdf,Jl0,Jl1,Jlog,Jof,Jhs,lam_s,recon] = computeObjMCDL(outputs);
        Jerr_true = sum(vec(abs(squeeze(recon)-squeeze(y_true)).^2))/2;
        obj_array(nn,r) = Jfn;
        true_error(nn,r) = Jerr_true;
        data_error(nn,r) = Jdf;
        l1_norm(nn,r) = Jl1;
        l0_norm(nn,r) = Jl0;
        log_penalty(nn,r) = Jlog;
    end
end

