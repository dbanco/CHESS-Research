function wrap_mmpad_norm2_FISTA_Circulant(datadir,P,outputdir)
%wrap_FISTA_Circulant 

%% load polar image
str2 = sprintf('%i',P.img);
fileName = ['mmpad_img_',str2,'.mat'];
fileDir = fullfile(datadir,fileName);
load(fileDir)

P.num_rad = size(polar_image,1);
P.num_theta = size(polar_image,2);

%% Generate unshifted basis function matrices
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);

%% Initialize solution
[m,n,t,r] = size(A0ft_stack);
x_init = ones(m,n,t,r);

%% call function
[x_hat, err, obj, l_0] = FISTA_Circulant(A0ft_stack,polar_image,x_init,P.params);


%% save output
save(fullfile(outputdir,sprintf('fista_fit_%i_%i.mat',P.ring_num,P.img)),'x_hat','err','polar_image','P')

end

