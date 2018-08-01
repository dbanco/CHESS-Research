%% Show mmpad video
fDir = 'D:\MMPAD_data\mmpad_ring1_fit';
for ring_num = 1
    for img_num = 540
        
        fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);
        load(fullfile(fDir,fName))
        polar_image = polar_image/norm(polar_image(:));
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
        img_fit = Ax_ft_2D(A0ft_stack,x_hat);
        lim1 = 0;
        lim2 = max(polar_image(:));
        % Plot both images
        figure(11)
        subplot(1,2,1)
        imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
        subplot(1,2,2)
        imshow(img_fit,'DisplayRange',[lim1 lim2],'Colormap',jet);
        title(sprintf('Error = %f',err(end)))
        pause(0.01)
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
for i = 1:P.num_var_t
    figure(i+10)
    for j = 1:P.num_var_r
        subplot(1,P.num_var_r,j)
        basis_func = shift2D(squeeze(A0_stack(:,:,i,j)),round(P.num_rad/2),round(P.num_theta/2));
        imshow(basis_func,'DisplayRange',[lim1 0.2],'Colormap',jet,'InitialMagnification','fit');
    end
    pause
end

%% View basis function
figure(2)
A0_stack = unshifted_basis_matrix_stack_norm2(P);
for i = 1:P.num_var_r
    figure(i+10)
    for j = 1:P.num_var_t
        subplot(1,P.num_var_t,j)
        basis_func = shift2D(squeeze(A0_stack(:,:,j,i)),round(P.num_rad/2),round(P.num_theta/2));
        imshow(basis_func,'DisplayRange',[lim1 0.2],'Colormap',jet,'InitialMagnification','fit');
    end
    pause
end

%% Compute time varying spread and error functions
fDir = 'D:\MMPAD_data\mmpad_ring1_fit';
az_spread = zeros(546,1);
rad_spread = zeros(546,1);
rel_err = zeros(546,1);
var_signal = zeros(546,P.num_var_t,P.num_var_r);
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
        rel_err(k) = err(end);
        var_signal_k = squeeze(sum(sum(x_hat,1),2));
        var_signal(k,:,:) =  var_signal_k;
        rad_var_signal = squeeze(sum(var_signal_k,2));
        az_var_signal = squeeze(sum(var_signal_k,1));
        var_sum = sum(var_signal_k(:));
        az_spread(k) = sqrt(P.var_rad)*az_var_signal'/var_sum;
        rad_spread(k) = sqrt(P.var_theta)*rad_var_signal/var_sum;
        
        k = k + 1;
    end
end

spreadDir = fullfile('D:','MMPAD_data','spread_results');
mkdir(spreadDir)
outFile = 'spread_mmpad_ring2_lambda_%2.2f.mat';
save(fullfile(spreadDir,sprintf(outFile,P.params.lambda)),'var_signal','rel_err','P')
        
%% Load spread data

load(fullfile(spreadDir,'spread_mmpad_ring1_lambda_0.10.mat'));
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

%% Show evolving var_signal
lim1 = 0;
lim2 = max(var_signal(:));
for i = 1:546
    imshow(squeeze(var_signal(i,:,:)),'DisplayRange',[lim1 lim2],'Colormap',jet,'InitialMagnification','fit');
    pause(0.2)
end