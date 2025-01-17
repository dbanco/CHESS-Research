function [true_error,obj_array,data_error,l1_norm,l0_norm,log_penalty,of_penalty,hs_penalty] = sim1_trials_results(sigmas,dataset,penalty,testType,topDir,num_trials)
%sim1_trials_results Summary of this function goes here
%   Detailed explanation goes here
NN = numel(sigmas)-1;
obj_array = zeros(NN,num_trials);
true_error = zeros(NN,num_trials);
data_error = zeros(NN,num_trials);
l1_norm = zeros(NN,num_trials);
l0_norm = zeros(NN,num_trials);
log_penalty = zeros(NN,num_trials);
of_penalty = zeros(NN,num_trials);
hs_penalty = zeros(NN,num_trials);

for nn = 2:9
    fprintf('Noise level %i \n',nn)
    [~,y_true,~,~,~] = gaus_example_switch_multiscale_dl(sigmas(nn),dataset);
    for r = 1:num_trials
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
        of_penalty(nn,r) = Jof;
        hs_penalty(nn,r) = Jhs;
    end
end

end