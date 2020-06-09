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
    parfor image_num = 1:num_ims
        xt_data = load(fullfile(input_dir,sprintf(baseFileName,1,image_num)));
        x_init = xt_data.x_hat;
        
        P_local = P;
        P_local.set = 1;
        P_local.params.lambda = lambda_values(image_num);
        P_local.params.time = image_num;
        
        % Reduce image to vector 
        try
            b = squeeze(sum(xt_data.polar_image,1));
        catch
            b = xt_data.polar_vector;
        end
        
        x_n = {};
        if image_num == 1
            xn_data = load(fullfile(input_dir,sprintf(baseFileName,1,2)),'x_hat')
            x_n{1} = 0;
            x_n{2} = xn_data.x_hat;
        elseif image_num == num_ims
            xn_data = load(fullfile(input_dir,sprintf(baseFileName,1,num_ims-1)),'x_hat')
            x_n{1} = xn_data.x_hat;
            x_n{2} = 0;
        else
            xn_data = load(fullfile(input_dir,sprintf(baseFileName,1,image_num-1)),'x_hat')
            x_n{1} = xn_data.x_hat;
            xn_data = load(fullfile(input_dir,sprintf(baseFileName,1,image_num+1)),'x_hat')
            x_n{2} = xn_data.x_hat;
        end
         
        % Coupled iterations
        [x_hat,err,obj] = convADMM_LASSO_Sherman_TVx_1D(A0ft_stack,b,x_init,x_n,P_local.params);  
        
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

