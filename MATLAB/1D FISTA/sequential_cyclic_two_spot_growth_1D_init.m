P.set = 1;
P.img = 1;

dataset = '/cluster/home/dbanco02/simulated_data_two_spot_growth_1D_25/';
output_dir = '/cluster/shared/dbanco02/two_spot_growth_1D_25_init1';
num_ims = 25;
mkdir(output_dir)
prefix = 'polar_image';

%% Universal Parameters
% Ring sampling parameters
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta= numel(polar_image);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'sqrt';
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;

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
params.gamma = 0.2;
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

for image_num = 1:num_ims
    im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
    %% Zero pad image
    b = zeroPad(im_data.polar_image,P.params.zeroPad);
    % Scale image by 2-norm
    b = b/norm(b(:));

    % Construct dictionary
    switch P.basis
        case 'norm2'
            A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
    end

    % Construct distance matrix
    N = P.num_var_t;
    THRESHOLD = 32;

    switch P.cost
            case 'l1'
                D = ones(N,N).*THRESHOLD;
                for i = 1:P.num_var_t
                    for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1]))
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

    x_init = zeros(size(A0ft_stack));
    for i = 1:P.num_var_t
            x_init(:,i) = b/P.num_var_t;
    end
  
    % Run FISTA updating solution and error array
    if image_num == 1
         [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D(A0ft_stack,b,x_init,P.params);
    else
         [x_hat, err, ~, ~,  obj, ~] = space_wasserstein_FISTA_Circulant_1D(A0ft_stack,b,vdfs,D,x_init,P.params);
    end
       
    new_vdf = squeeze(sum(x_hat,1))/sum(x_hat(:));
    vdfs = {new_vdf};
    
    save_output(output_dir,baseFileName,x_hat,err,im_data.polar_image,P,image_num);
    save_obj(output_dir,0,image_num,obj);
end

function save_output(output_dir,baseFileName,x_hat,err,polar_image,P,image_num)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','polar_image','P');
end
function save_obj(output_dir,pass,image_num,obj)
    save(fullfile(output_dir,sprintf('objective_%i_%i.mat',pass,image_num)),'obj');
end
