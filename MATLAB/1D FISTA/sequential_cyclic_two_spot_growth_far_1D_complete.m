P.set = 1;

dataset = '/cluster/home/dbanco02/simulated_data_two_spot_growth_1D_far';
num_ims = 25;
rescale = 100;

%% Universal Parameters
% Ring sampling parameters
prefix = 'polar_image';
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 15;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

%% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1000;
params.t_k = 1;
params.lambda = 0.03;
params.wLam = 25;
params.gamma = 0.075;
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


iii = 1;


init_dir = ['/cluster/shared/dbanco02/two_spot_growth_far_1D_init' num2str(iii)];
output_dirA = ['/cluster/shared/dbanco02/two_spot_growth_far_1D_' num2str(iii) 'a'];
output_dirB = ['/cluster/shared/dbanco02/two_spot_growth_far_1D_' num2str(iii) 'b'];

% init_dir = ['D:\CHESS_data\two_spot_growth_1D_25_renorm_init' num2str(iii)];
% output_dirA = ['D:\CHESS_data\two_spot_growth_1D_25_renorm_' num2str(iii) 'a'];
% output_dirB = ['D:\CHESS_data\two_spot_growth_1D_25_renorm_' num2str(iii) 'b'];

mkdir(init_dir)
mkdir(output_dirA)
mkdir(output_dirB)

for jjj = 1
    % setup io directories
    if jjj == 1
        input_dir = init_dir;
        output_dir = init_dir;
    elseif jjj == 2
        input_dir = init_dir;
        output_dir = output_dirA;
    elseif mod(jjj,2)
        input_dir = output_dirA;
        output_dir = output_dirB;
    else
        input_dir = output_dirB;
        output_dir = output_dirA;
    end
    % extract vdfs if beyond first pass
    if jjj > 1
        vdf_array = cell(num_ims,1);
        for ii = 1:num_ims
            f_data = load(fullfile(input_dir,sprintf(baseFileName,1,ii)));
            vdf_array{ii} = squeeze(sum(f_data.x_hat,1))/sum(f_data.x_hat(:));
        end
        new_vdf_array = cell(num_ims,1);
    end
    vdfs = {};
    % iterate over each image
    for image_num = 1:num_ims
        im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
        %% Zero pad image
        b = zeroPad(im_data.polar_image,P.params.zeroPad);
        % Scale image by 2-norm
        b = b/rescale;
        
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
        
        x_init = zeros(size(A0ft_stack));
        for i = 1:P.num_var_t
            x_init(:,i) = b/P.num_var_t;
        end
        
        if jjj == 1
            % First iteration initializes causally
            if image_num == 1
                [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D(A0ft_stack,b,x_init,P.params);
            else
                [x_hat, err, ~, ~,  obj, ~] = space_wasserstein_FISTA_Circulant_1D(A0ft_stack,b,vdfs,D,x_init,P.params);
            end
            new_vdf = squeeze(sum(x_hat,1))/sum(x_hat(:));
            vdfs = {new_vdf};
            
        else
            % Consecutive iterations use t-1,t+1
            if image_num == 1
                vdfs = vdf_array(2);
            elseif image_num == num_ims
                vdfs = vdf_array(num_ims-1);
            else
                vdfs = {vdf_array{image_num-1},vdf_array{image_num+1}};
            end
            [x_hat, err, ~, ~,  obj, ~] = space_wasserstein_FISTA_Circulant_1D(A0ft_stack,b,vdfs,D,x_init,P.params);
            
            new_vdf_array{image_num} = squeeze(sum(x_hat,1))/sum(x_hat(:));
        end
        
        % Output data
        save_output(output_dir,baseFileName,x_hat,err,im_data.polar_image,P,image_num);
        save_obj(output_dir,jjj,image_num,obj);
    end
    if jjj > 1
        vdf_array = new_vdf_array;
    end
end
parpool(13)
parfor jjj = 2:11
    % setup io directories
    if jjj == 1
        input_dir = init_dir;
        output_dir = init_dir;
    elseif jjj == 2
        input_dir = init_dir;
        output_dir = output_dirA;
    elseif mod(jjj,2)
        input_dir = output_dirA;
        output_dir = output_dirB;
    else
        input_dir = output_dirB;
        output_dir = output_dirA;
    end
    % extract vdfs if beyond first pass
    if jjj > 1
        vdf_array = cell(num_ims,1);
        for ii = 1:num_ims
            f_data = load(fullfile(input_dir,sprintf(baseFileName,1,ii)));
            vdf_array{ii} = squeeze(sum(f_data.x_hat,1))/sum(f_data.x_hat(:));
        end
        new_vdf_array = cell(num_ims,1);
    end
    vdfs = {};
    % iterate over each image
    for image_num = 1:num_ims
        im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
        %% Zero pad image
        b = zeroPad(im_data.polar_image,P.params.zeroPad);
        % Scale image by 2-norm
        b(b<0) = 0;
        b = b/rescale;
        
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
        
        x_init = zeros(size(A0ft_stack));
        for i = 1:P.num_var_t
            x_init(:,i) = b/P.num_var_t;
        end
        
        if jjj == 1
            % First iteration initializes causally
            if image_num == 1
                [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D(A0ft_stack,b,x_init,P.params);
            else
                [x_hat, err, ~, ~,  obj, ~] = space_wasserstein_FISTA_Circulant_1D(A0ft_stack,b,vdfs,D,x_init,P.params);
            end
            new_vdf = squeeze(sum(x_hat,1))/sum(x_hat(:));
            vdfs = {new_vdf};
            
        else
            % Consecutive iterations use t-1,t+1
            if image_num == 1
                vdfs = vdf_array(2);
            elseif image_num == num_ims
                vdfs = vdf_array(num_ims-1);
            else
                vdfs = {vdf_array{image_num-1},vdf_array{image_num+1}};
            end
            [x_hat, err, ~, ~,  obj, ~] = space_wasserstein_FISTA_Circulant_1D(A0ft_stack,b,vdfs,D,x_init,P.params);
            
            new_vdf_array{image_num} = squeeze(sum(x_hat,1))/sum(x_hat(:));
        end
        
        % Output data
        save_output(output_dir,baseFileName,x_hat,err,im_data.polar_image,P,image_num);
        save_obj(output_dir,jjj,image_num,obj);
    end
    if jjj > 1
        vdf_array = new_vdf_array;
    end
end


function save_output(output_dir,baseFileName,x_hat,err,polar_image,P,image_num)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','polar_image','P');
end
function save_obj(output_dir,pass,image_num,obj)
    save(fullfile(output_dir,sprintf('objective_%i_%i.mat',pass,image_num)),'obj');
end
