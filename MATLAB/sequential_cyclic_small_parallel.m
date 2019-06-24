P.set = 1;
P.img = 1;

% dataset = 'D:\CHESS_data\simulated_data_small';
% output_dir = 'D:\CHESS_data\small_wass_test_results';
dataset = '/cluster/home/dbanco02/simulated_data_two_phase/';
output_dir = '/cluster/shared/dbanco02/seq_two_phase3';
mkdir(output_dir)
prefix = 'polar_image';

%% Universal Parameters
% Ring sampling parameters
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [20,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

%% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1000;
params.t_k = 1;
params.lambda = 0.0359;
params.wLam = 25;
params.gamma = 1;
params.beta = 1.2;
params.maxIter = 800;
params.maxIterReg = 800;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

baseFileName = 'fista_fit_%i_%i.mat';

for jjj = 1:20
    parfor image_num = 1:20
        f_data = load(fullfile(output_dir,sprintf(baseFileName,1,image_num));
        %% Zero pad image
        im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
        P = f_data.P
        b = zeroPad(m_data.polar_image,P.params.zeroPad);
        % Scale image by 2-norm
        b = b/norm(b(:));
        P.num_rad = size(b,1);
        P.num_theta = size(b,2);

        % Construct dictionary
        switch f_data.P.basis
            case 'norm2'
                A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
            case 'norm1'
                A0ft_stack = unshifted_basis_matrix_ft_stack_norm(P);
            case 'max'
                A0ft_stack = unshifted_basis_matrix_ft_stack(P);
        end

        % Construct distance matrix
        N = P.num_var_t*P.num_var_r;
        THRESHOLD = 32;

        D = ones(N,N).*THRESHOLD;
        for i = 1:P.num_var_t
            for j = 1:P.num_var_r
                for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1])
                    for jj = max([1 j-THRESHOLD+1]):min([P.num_var_r j+THRESHOLD-1])
                        ind1 = i + (j-1)*P.num_var_t;
                        ind2 = ii + (jj-1)*P.num_var_t;
                        D(ind1,ind2)= sqrt((i-ii)^2+(j-jj)^2); 
                    end
                end
            end
        end
        D = D./max(D(:));

        x_init = f_data.x_hat

        % Run FISTA updating solution and error array
        vdfs = load_neighbors_vdf(output_dir,baseFileName,P);
        [x_hat, err, ~, ~,  ~, ~] = space_wasserstein_FISTA_Circulant(A0ft_stack,b,vdfs,D,x_init,P.params);

        save_output(output_dir,baseFileName,x_hat,err,polar_image,P);
    end
end
function save_output(output_dir,baseFileName,x_hat,err,polar_image,P)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,P.img)),'x_hat','err','polar_image','P');
end

% % Compare error/ wasserstein distance/ vdfs
% figure(11)
% subplot(1,3,1)
% imagesc(vdf_true)
% [awmv_az_true,awmv_rad_true] = computeAWMV_vdf(vdf_true,P.var_theta,P.var_rad);
% title(['AWMV \theta: ',num2str(awmv_az_true),'   AWMV r: ',num2str(awmv_rad_true)])
% 
% subplot(1,3,2)
% vdf_unreg = squeeze(sum(sum(x_unreg,1),2))/sum(x_unreg(:));
% Wd_unreg = sinkhornKnoppTransport(vdf_true, vdf_unreg, P.params.wLam, D);
% [awmv_az_unreg,awmv_rad_unreg] = computeAWMV_vdf(vdf_unreg,P.var_theta,P.var_rad);
% imagesc(vdf_unreg)
% 
% title(['Unreg Wd: ',num2str(Wd_unreg),'   AWMV \theta: ',num2str(awmv_az_unreg),'   AWMV r: ',num2str(awmv_rad_unreg)])
% 
% subplot(1,3,3)
% vdf_reg = squeeze(sum(sum(x_hat,1),2))/sum(x_hat(:));
% Wd_reg = sinkhornKnoppTransport(vdf_true, vdf_reg, P.params.wLam, D);
% [awmv_az_reg,awmv_rad_reg] = computeAWMV_vdf(vdf_reg,P.var_theta,P.var_rad);
% imagesc(vdf_reg)
% title(['Reg Wd: ',num2str(Wd_reg),'   AWMV \theta: ',num2str(awmv_az_reg),'   AWMV r: ',num2str(awmv_rad_reg)])
% 
% figure(22)
% imagesc(polar_image)
% 
% 
% 
