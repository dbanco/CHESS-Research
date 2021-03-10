function wrap_convADMM_LASSO_ISM_Mask_1D(dataset,P,output_dir)
%wrap_convADMM_LASSO_CG_TVphi_1D Wrapper function to call 
%   convADMM_LASSO_CG_TVphi_1D during in parallel processing

N = P.num_theta;
K = P.num_var_t;
T = P.num_ims;

% Construct dictionary
P.num_theta = N + 2*zPad;
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
P.num_theta = N;
% Load data
B = zeros(N,T);
for j = 1:T
  b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
    % Reduce image to vector if needed
    try
        b = P.dataScale*sum(b_data.polar_image,1);
        b(129:133,j) = (b(128) + b(134))/2;
    catch
        b = P.dataScale*b_data.polar_vector(1:179);
    end
    b = zeroPad(b,zPad);
    B(:,j) = b;

end

% Init solution
x_init = zeros(N,K);
for t = 1:T
    % Solve
    [x_hat,err,obj,l1_norm] = convADMM_LASSO_MaskDecouple_1D(A0ft_stack,B(:,t),x_init,P.params); 

    % Output data
    save(fullfile(output_dir,sprintf(P.baseFileName,P.set,t)),...
        'x_hat','err','obj','l1_norm','P');
end

end

