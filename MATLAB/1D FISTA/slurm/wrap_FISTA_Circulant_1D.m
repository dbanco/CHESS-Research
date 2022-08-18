function wrap_FISTA_Circulant_1D(datadir,P,outputdir)
%wrap_norm2_FISTA_Circulant 

%% load polar image
try
    str1 = sprintf('%i',P.set);
    str2 = sprintf('%i',P.img);
    fileDir = [datadir,'_',str1,'_',str2,'.mat'];
    load(fileDir)
catch
    str1 = sprintf('%i',P.img);
    fileDir = [datadir,'_',str1,'.mat'];
    load(fileDir)
end

%% Zero pad image
b = zeroPad(polar_vector,P.params.zeroPad);

P.num_theta = numel(b);

%% Generate unshifted basis function matrices
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
    case 'norm1'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm(P);
    case 'max'
        A0ft_stack = unshifted_basis_vector_ft_stack(P);
end

% Initialize solution
x_init = rand(size(A0ft_stack));
x_init = x_init/norm(x_init(:));
x_init = forceMaskToZeroArray(x_init,P.params.zeroMask);

% Scale image by 2-norm
b = b/norm(b(:));

%% Call function
[x_hat, err, obj, l_0] = FISTA_Circulant_1D(A0ft_stack,b,x_init,P.params);

%% Save output
save(fullfile(outputdir,sprintf('fista_fit_%i_%i.mat',P.set,P.img)),'x_hat','err','polar_vector','P')

end

