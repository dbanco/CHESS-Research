clear all
close all

baseDir = 'D:\CHESS_data\';

datasetArray = {'simulated_two_spot_1D_noise2_independent_8',...
                'simulated_two_spot_1D_noise2_independent_quadratic_8',...
                'noise2_1D_8_1b',...
                'noise_1D_Quadratic_8_1b'};

num_imgs = 10;

figure(222)
[ha1, pos1] = tight_subplot(4,10,[.005 .005],[.01 .01],[.01 .01]); 
figure(223)
[ha2, pos2] = tight_subplot(4,10,[.1 .02],[.05 .05],[.05 .05]); 
kk = 1;

for i = 1:numel(datasetArray)
    
    fDir = [baseDir,datasetArray{i}];
    
    for img_num = 1:num_imgs
        fprintf('Image %i\n',img_num)
        fName = sprintf('fista_fit_%i_%i.mat',1,img_num);
        load(fullfile(fDir,fName))
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
        img_fit = Ax_ft_1D(A0ft_stack,x_hat);
        vdf_az = squeeze(sum(x_hat,1))./sum(x_hat(:));
        try
            axes(ha1(kk));
            
            % Plot fits
            hold on
            plot(polar_image./norm(polar_image(:)))
            plot(img_fit)
            
            axes(ha2(kk));
            % Plot vdfs
            plot(vdf_az)
            title(sprintf('awmv = %0.2f',sqrt(P.var_theta)*vdf_az'))
            kk = kk + 1;
            
        catch
                
        end
        
    %     if k > 2
    %         wass_dist(k) = sinkhornKnoppTransport(var_signal_k(:),vdf_last(:),P.params.wLam,D);
    %     end
    end
end