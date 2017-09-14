function wrap_FISTA_Circulant(datadir,P,outputdir)
%wrap_FISTA_Circulant 

%% Generate unshifted basis function matrices
if P.weight
    A0ft_stack = unshifted_basis_matrix_ft_stack_weight(P.var_theta,P.var_rad,P.dtheta,P.drad,P.num_theta,P.num_rad,P.betap);
else
    A0ft_stack = unshifted_basis_matrix_ft_stack(P.var_theta,P.var_rad,P.dtheta,P.drad,P.num_theta,P.num_rad);
end

%% load polar image
str1 = sprintf('%i',P.load_step);
str2 = sprintf('%i',P.img);
fileName = ['polar_image_',str1,'_',str2,'.mat'];
fileDir = fullfile(datadir,fileName);
load(fileDir)

%% call function
[x_hat, err, obj, l_0] = FISTA_Circulant(A0ft_stack,polar_image,P.params);

%% save output
save(fullfile(outputdir,sprintf('fista_fit_%i_%i_betap_%2.2f.mat',P.load_step,P.img,P.betap)),'x_hat','err','polar_image','P')
end
