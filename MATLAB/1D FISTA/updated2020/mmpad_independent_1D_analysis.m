clear all
close all

baseDir = 'D:\CHESS_data\';
baseDataset = 'mmpad_1D_subset_independent';
suffix = '';
fitName = '';

gamma_vals = [0.1]; %[0.005 0.01 0.02 0.05 0.1]; % 0.2 0.5 1];
num_imgs = 10;

fDir = [baseDir,baseDataset,suffix];
fName = sprintf('fista_fit_%i_%i.mat',1,1);
load(fullfile(fDir,fName))

az_spread = zeros(num_imgs,numel(gamma_vals));
rel_err = zeros(num_imgs,numel(gamma_vals));
sparsity = zeros(num_imgs,numel(gamma_vals));
objective = zeros(num_imgs,numel(gamma_vals));
var_signal = zeros(num_imgs,numel(gamma_vals),P.num_var_t);
prefix = 'polar_image';


for gam_num = 1
    
    fDir = [baseDir,baseDataset,suffix];

    for img_num = 1:num_imgs
        k = img_num;
        fprintf('Image %i\n',k)
        fName = sprintf('fista_fit_%i_%i.mat',1,img_num);
        load(fullfile(fDir,fName))
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
        img_fit = Ax_ft_1D(A0ft_stack,x_hat);
        rel_err(k,gam_num) = err(end-1);
        var_signal_k = squeeze(sum(x_hat,1));
        var_signal(k,gam_num,:) =  var_signal_k;
        az_var_signal = squeeze(sum(var_signal_k,1));
        var_sum = sum(var_signal_k(:));
        az_spread(k,gam_num) = sum(sqrt(P.var_theta(:)).*az_var_signal(:))/var_sum;
        sparsity(k,gam_num) = sum(x_hat(:)>0);
        objName = sprintf('objective_%i_%i.mat',0,img_num);
        load(fullfile(fDir,objName))
        objective(k,gam_num) = obj(end-1);
    %     if k > 2
    %         wass_dist(k) = sinkhornKnoppTransport(var_signal_k(:),vdf_last(:),P.params.wLam,D);
    %     end

        k = k + 1;
        vdf_last = var_signal_k;
    end
end

spreadDir = fullfile('D:','CHESS_data','spread_results');
outFile = ['spread_',baseDataset,fitName,'.mat'];
save(fullfile(spreadDir,outFile),...
    'var_signal','rel_err','P','az_spread','rel_err','sparsity','objective')%,'wass_dist')
 