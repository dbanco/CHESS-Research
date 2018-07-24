%% Show mmpad video
fDir = 'D:\MMPAD_data\mmpad_ring1_fit1';
for ring_num = 1
    for img_num = 1:546
        
        fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);
        load(fullfile(fDir,fName))
        polar_image = polar_image/norm(polar_image(:));
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
        img_fit = Ax_ft_2D(A0ft_stack,x_hat);
        lim1 = 0;
        lim2 = max(polar_image(:));
        % Plot both images
        figure(1)
        subplot(1,2,1)
        imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
        subplot(1,2,2)
        imshow(img_fit,'DisplayRange',[lim1 lim2],'Colormap',jet);
        pause
    end
end

%% Load mmpad fit
img_num = 25;
ring_num = 1;

fDir = 'D:\MMPAD_data\mmpad_ring1_fit4';
fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);

load(fullfile(fDir,fName))
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
img_fit = Ax_ft_2D(A0ft_stack,x_hat);

lim1 = 0;
lim2 = 300000;
% Plot both images
figure(1)
subplot(1,2,1)
imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
subplot(1,2,2)
imshow(img_fit,'DisplayRange',[lim1 lim2],'Colormap',jet);

%% View basis function
figure(2)
A0_stack = unshifted_basis_matrix_stack_norm2(P);
for i = 1:10
    for j = 1:10
        subplot(10,10,(i-1)*10+j)
        basis_func = shift2D(squeeze(A0_stack(:,:,i,j)),round(P.num_rad/2),round(P.num_theta/2));
        imshow(basis_func','DisplayRange',[lim1 0.2],'Colormap',jet);
    end
end

%% Compute time varying spread and error functions
fDir = 'D:\MMPAD_data\mmpad_ring1_fit';
az_spread = zeros(546,1);
rad_spread = zeros(546,1);
rel_err = zeros(546,1);
k = 1;
for ring_num = 1
    for img_num = 1:546
        fprintf('Image %i\n',k)
        fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);
        load(fullfile(fDir,fName))
%         P.num_rad = size(polar_image,1);   WRONG
%         P.num_theta = size(polar_image,2); DIDNT WORK
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
        img_fit = Ax_ft_2D(A0ft_stack,x_hat);
        rel_err(k) = norm(polar_image(:)-img_fit(:))/norm(polar_image(:));
        var_signal = squeeze(sum(sum(x_hat,1),2));
        rad_var_signal = squeeze(sum(var_signal,2));
        az_var_signal = squeeze(sum(var_signal,1));
        var_sum = sum(var_signal(:));
        az_spread(k) = sqrt(P.var_rad)*az_var_signal'/var_sum;
        rad_spread(k) = sqrt(P.var_theta)*rad_var_signal/var_sum;
        
        k = k + 1;
    end
end

%% Plot time varying spread and error functions
figure(1)
plot(az_spread,'-')
hold on
plot(rad_spread,'-')
legend('azimuthal','radial','Location','Best')
title('Spread Metrics')
xlabel('Time')
ylabel('MVar')

figure(2)
plot(rel_err,'-')
title('Relative Error')
xlabel('Time')
