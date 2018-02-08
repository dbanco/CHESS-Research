function init_FISTA_Circulant( data_dir,P,output_dir )
%init_FISTA_Circuleant Runs FISTA_Circulant loading input files and saving
% ouput files

% Load image file
str1 = sprintf('%i',P.load_step);
str2 = sprintf('%i',P.img);
fileName = ['polar_image_',str1,'_',str2,'.mat'];
fileDir = fullfile(data_dir,fileName);
load(fileDir)

% Construct dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack(P);

% Initialize solution
[m,n,t,r] = size(A0ft_stack);
x_init = ones(m,n,t,r);

%% Run FISTA
[x_hat, err, obj, l_0] = FISTA_Circulant(A0ft_stack,polar_image,x_init,P.params);

%% Save outputs
save(fullfile(output_dir,sprintf('spatial_fit_%i_%i.mat',P.load_step,P.img,P.task)),'x_hat','err','polar_image','P')


end

