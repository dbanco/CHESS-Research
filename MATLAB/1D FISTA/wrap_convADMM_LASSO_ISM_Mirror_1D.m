function wrap_convADMM_LASSO_ISM_Mirror_1D(dataset,P,output_dir)
%wrap_convADMM_LASSO_CG_TVphi_1D Wrapper function to call 
%   convADMM_LASSO_CG_TVphi_1D during in parallel processing

N = P.num_theta;
K = P.num_var_t;
T = P.num_ims;

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Load data
B = zeros(N,T);
for j = 1:T
  b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
    % Reduce image to vector if needed
    try
        b = sum(b_data.polar_image,1);
        b(129:133) = (b(128) + b(134))/2;
    catch
        b = b_data.polar_vector;
    end

    % Mirror data
    nn = numel(b);
    pad1 = floor(nn/2);
    pad2 = ceil(nn/2);
    N = nn + pad1 + pad2;
    b_mirror = zeros(N,1);
    b_mirror((pad1+1):(pad1+nn)) = b;
    b_mirror((1+N-pad2):N) = flip(b((nn-pad2+1):nn));
    b_mirror(1:pad1) = flip(b(1:pad1));    
    B(:,j) = b_mirror;

end

% Rescale data
P.dataScale = 1/mean(B(:));
B = B*P.dataScale;

% Init solution
x_init = zeros(N,K);
for t = 1:T
    % Solve
    bnorm = norm(B(:,t));
    [x_hat,err,obj,l1_norm,rho] = convADMM_LASSO_Sherman_1D(A0ft_stack/bnorm,B(:,t)/bnorm,x_init,P.params); 

    % Output data
    save(fullfile(output_dir,sprintf(P.baseFileName,P.set,t)),...
        'x_hat','err','obj','l1_norm','rho','P');
end

end

