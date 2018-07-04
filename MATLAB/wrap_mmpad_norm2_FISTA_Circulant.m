function wrap_mmpad_norm2_FISTA_Circulant(datadir,P,outputdir)
%wrap_FISTA_Circulant 

%% Generate unshifted basis function matrices
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);

%% load polar image
str2 = sprintf('%i',P.img);
fileName = ['mmpad_img_',str2,'.mat'];
fileDir = fullfile(datadir,fileName);
load(fileDir)

% Initialize solution
[m,n,t,r] = size(A0ft_stack);
x_init = ones(m,n,t,r);

%% call function
[x_hat, err, obj, l_0] = FISTA_Circulant(A0ft_stack,polar_image,x_init,P.params);


%% save output
save(fullfile(outputdir,sprintf('fista_fit_%i_%i.mat',P.ring_num,P.img)),'x_hat','err','polar_image','P')

end

