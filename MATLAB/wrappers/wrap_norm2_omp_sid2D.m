function wrap_norm2_omp_sid2D(datadir,P,outputdir)
%wrap_norm2_omp_sid2D 

%% load polar image
str1 = sprintf('%i',P.load_step);
str2 = sprintf('%i',P.img);
fileName = ['polar_image_',str1,'_',str2,'.mat'];
fileDir = fullfile(datadir,fileName);
load(fileDir)


%% call function
x_hat = saSID2Domp(polar_image, P.sid, P.s,'mexOMP',P.bSize,0);
b_fit = multSID2D(P.sid, x_hat, P.num_rad, P.num_theta);
err = norm(polar_image-b_fit)/norm(polar_image);

%% save output
save(fullfile(outputdir,sprintf('fista_fit_%i_%i.mat',P.load_step,P.img)),'x_hat','err','polar_image','P')

end

