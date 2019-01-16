function wrap_GCMP(datadir,P,outputdir)
%wrap_GCMP

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
b = zeroPad(polar_image,P.params.zeroPad);

P.num_rad = size(b,1);
P.num_theta = size(b,2);

%% Generate unshifted basis function matrices
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
    case 'norm1'
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm(P);
    case 'max'
        A0ft_stack = unshifted_basis_matrix_ft_stack(P);
end

%% Call function
x_hat = GCMP(A0ft_stack,b,P.params);

%% Save output
save(fullfile(outputdir,sprintf('gcmp_fit_%i_%i.mat',P.set,P.img)),'x_hat','polar_image','P')

end

