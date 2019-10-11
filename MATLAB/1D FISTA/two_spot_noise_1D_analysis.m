clear all
close all

baseDir = 'D:\CHESS_data\';
baseDataset = 'simulated_two_spot_1D_noise_independent_quadratic';

num_imgs = 10;

fDir = [baseDir,baseDataset,'_1'];
fName = sprintf('fista_fit_%i_%i.mat',1,1);
load(fullfile(fDir,fName))

az_spread = zeros(num_imgs,12);
rel_err = zeros(num_imgs,12);
sparsity = zeros(num_imgs,12);
objective = zeros(num_imgs,12);
var_signal = zeros(num_imgs,12,P.num_var_t);
prefix = 'polar_image';


% Construct distance matrix
N = P.num_var_t;
THRESHOLD = 32;

switch P.cost
    case 'l1'
        D = ones(N,N).*THRESHOLD;
        for i = 1:P.num_var_t
            for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1])
                D(i,ii)= abs(i-ii); 
            end
        end

        D = D./max(D(:));

    case 'l2'
        D = ones(N,N).*THRESHOLD;
        for i = 1:P.num_var_t
            for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1])
                D(i,ii)= (i-ii)^2; 
            end
        end

        D = D./max(D(:));

    case 'wass'
        D = ones(N,N).*THRESHOLD;
        for i = 1:P.num_var_t
            for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1])                        
                D(i,ii)= P.var_theta(i) + P.var_theta(ii) -...
                         2*sqrt(P.var_theta(i)*P.var_theta(ii));
            end
        end
        D = D./max(D(:));

    case 'sqrt'
        D = ones(N,N).*THRESHOLD;
        for i = 1:P.num_var_t
            for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1])
                D(i,ii)= sqrt(abs(i-ii));
            end
        end
        D = D./max(D(:));
end

for j = 1:12
    
    suffix = sprintf('_%i',j);
    
    fDir = [baseDir,baseDataset,suffix];

    for img_num = 1:num_imgs
        k = img_num;
        fprintf('Image %i\n',k)
        fName = sprintf('fista_fit_%i_%i.mat',1,img_num);
        load(fullfile(fDir,fName))
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
        img_fit = Ax_ft_1D(A0ft_stack,x_hat);
        rel_err(k,j) = err(end-1);
        var_signal_k = squeeze(sum(x_hat,1));
        var_signal(k,j,:) =  var_signal_k;
        az_var_signal = squeeze(sum(var_signal_k,1));
        var_sum = sum(var_signal_k(:));
        az_spread(k,j) = sum(sqrt(P.var_theta(:)).*az_var_signal(:))/var_sum;
        sparsity(k,j) = sum(x_hat(:)>0);
        objName = sprintf('objective_%i_%i.mat',0,img_num);
        load(fullfile(fDir,objName))
        objective(k,j) = obj(end-1);
    %     if k > 2
    %         wass_dist(k) = sinkhornKnoppTransport(var_signal_k(:),vdf_last(:),P.params.wLam,D);
    %     end

        k = k + 1;
        vdf_last = var_signal_k;
    end
end

spreadDir = fullfile('D:','CHESS_data','spread_results');
outFile = ['spread_',baseDataset,'.mat'];
save(fullfile(spreadDir,outFile),...
    'var_signal','rel_err','P','az_spread','rel_err','sparsity','objective')%,'wass_dist')
 