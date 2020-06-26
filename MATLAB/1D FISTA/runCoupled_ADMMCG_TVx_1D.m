function runCoupled_ADMMCG_TVx_1D( P, Pc )
%runCoupledISM_TVx_1D

% Unpack parameters
P.params.rho1 = Pc.rho1;
P.params.lambda1 = Pc.lambda1;
P.params.rho2 = Pc.rho2;
P.params.lambda2 = Pc.lambda2;
P.params.maxIter = Pc.maxIterReg;
P.params.tolerance = Pc.tolerance;

dataset = Pc.dataset;
% lambda_values = Pc.lambda_values;
num_outer_iters = Pc.num_outer_iters;
init_dir = Pc.init_dir;
output_dirA = Pc.output_dirA;
output_dirB = Pc.output_dirB;
output_dirFinal = Pc.output_dirFinal;
baseFileName = Pc.baseFileName;
num_ims = Pc.num_ims;
prefix = Pc.prefix;

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
    
for jjj = 1:num_outer_iters
    % Setup io directories
    if jjj == 1
        input_dir = init_dir;
        output_dir = output_dirA;
    elseif ~mod(jjj,2)
        input_dir = output_dirA;
        output_dir = output_dirB;
    else
        input_dir = output_dirB;
        output_dir = output_dirA;
    end
    

    % Load independent solution
    xt_data = load(fullfile(input_dir,sprintf(baseFileName,1)));
    X_init = xt_data.X_hat;

    % Track image number and lambda parameter
    P_local = P;
    P_local.set = 1;
%     P_local.params.lambda2 = lambda_values(image_num);

    B = zeros(P.num_theta,num_ims);
    for i = 1:num_ims
      b_data = load(fullfile(dataset,[prefix,'_',num2str(i),'.mat']));
        % Reduce image to vector if needed
        try
            b = squeeze(sum(b_data.polar_image,1));
        catch
            b = b_data.polar_vector;
        end
        B(:,i) = b';
    end

    % Solve
    [X_hat,err,obj,~,~] = convADMM_LASSO_CG_TVx_1D(A0ft_stack,B,X_init,P_local.params);  

    % Output data
    try
        save_output(output_dir,X_hat,err,P_local,Pc);
        save_obj(output_dir,jjj,obj);
    catch
        save_output(output_dir,X_hat,err,P_local,Pc);
        save_obj(output_dir,jjj,obj);
    end
end

% Move solution to final output dir
movefile(fullfile(output_dir,sprintf(baseFileName,P.set)),output_dirFinal);

end

function save_output(output_dir,X_hat,err,P,Pc)
    save(fullfile(output_dir,sprintf(Pc.baseFileName,P.set)),'X_hat','err','P','Pc');
end
function save_obj(output_dir,pass,obj)
    save(fullfile(output_dir,sprintf('objective_%i_%i.mat',pass)),'obj');
end

