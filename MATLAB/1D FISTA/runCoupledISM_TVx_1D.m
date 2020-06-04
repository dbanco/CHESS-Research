function runCoupledISM_TVx_1D( P, Pc )
%runCoupledISM_TVx_1D

% Unpack parameters
P.params.rho2 = Pc.rho2;
P.params.lambda2 = Pc.lambda2;
P.params.maxIterReg = Pc.maxIterReg;

dataset = Pc.dataset;
lambda_values = Pc.lambda_values;
num_outer_iters = Pc.num_outer_iters;
init_dir = Pc.init_dir;
output_dirA = Pc.output_dirA;
output_dirB = Pc.output_dirB;
baseFileName = Pc.baseFileName;
num_ims = Pc.num_ims;
prefix = Pc.prefix;

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end

if Pc.preInitialized
   start_ind = Pc.preInitialized;
else
    start_ind = 1;
end
    
for jjj = start_ind:num_outer_iters
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
    
    % iterate over each image
    for image_num = 1:num_ims
        im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
        n_ind = [image_num-1, image_num+1];
        n_ind = n_ind( (n_ind>=1)&(n_ind<=num_ims) );
        x_n = {};
        k = 1;
        for i = n_ind
            x_data = load(fullfile(input_dir,sprintf(baseFileName,1,i)),'x_hat')
            x_n{k} = x_data.x_hat;
            k = k + 1;
        end

        % Reduce image to vector 
        try
            b = squeeze(sum(im_data.polar_image,1));
        catch
            b = im_data.polar_vector;
        end
        
        P_local = P;
        P_local.set = 1;
        P_local.params.plotProgress = 1;
        % Use selected lambda
        P_local.params.lambda = lambda_values(image_num);
 
        x_init = zeros(size(A0ft_stack));
        for i = 1:P_local.num_var_t
            x_init(:,i) = zeroPad(b/P_local.num_var_t,P_local.params.zeroPad);
        end
        
        if jjj == 1
            % Initialization 
            switch Pc.initialization
                case 'causal'
                    if image_num == 1
                        [x_hat,err,obj] = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_init,P_local.params);
                    else
                        [x_hat,err,obj] = convADMM_LASSO_Sherman_TVx_1D(A0ft_stack,b,x_init,x_n,P_local.params);
                    end
                case 'simultaneous'
                    [x_hat,err,obj] = convADMM_LASSO_Sherman_1D(A0ft_stack,bn,x_init,P_local.params);
            end     
        else
            [x_hat,err,obj] = convADMM_LASSO_Sherman_TVx_1D(A0ft_stack,b,x_init,x_n,P_local.params);  
        end
        
        % Output data
        try
            save_output(output_dir,baseFileName,x_hat,err,im_data.polar_image,P_local,Pc,image_num);
            save_obj(output_dir,jjj,image_num,obj);
        catch
            save_output(output_dir,baseFileName,x_hat,err,im_data.polar_vector,P_local,Pc,image_num);
            save_obj(output_dir,jjj,image_num,obj);
        end
    end
end

end

function save_output(output_dir,baseFileName,x_hat,err,polar_image,P,Pc,image_num)
    save(fullfile(output_dir,sprintf(baseFileName,P.set,image_num)),'x_hat','err','polar_image','P','Pc');
end
function save_obj(output_dir,pass,image_num,obj)
    save(fullfile(output_dir,sprintf('objective_%i_%i.mat',pass,image_num)),'obj');
end

