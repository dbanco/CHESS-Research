%% Compute time varying spread and error functions
fDir = 'D:\CHESS_data\synth2_norm2_correct1';
fName = sprintf('fista_fit_%i_%i.mat',0,0);
load(fullfile(fDir,fName))
az_spread = zeros(5,5);
rad_spread = zeros(5,5);
rel_err = zeros(5,5);
var_signal = zeros(5,5,P.num_var_t,P.num_var_r);
k = 0;
for row = 1:5
    for col = 1:5
        fprintf('Image %i\n',k)
        fName = sprintf('fista_fit_%i_%i.mat',0,k);
        load(fullfile(fDir,fName))
%         P.num_rad = size(polar_image,1);   WRONG
%         P.num_theta = size(polar_image,2); DIDNT WORK
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
        img_fit = Ax_ft_2D(A0ft_stack,x_hat);
        rel_err(row,col) = err(end);
        var_signal_k = squeeze(sum(sum(x_hat,1),2));
        var_signal(row,col,:,:) =  var_signal_k;
        az_var_signal = squeeze(sum(var_signal_k,2));
        rad_var_signal = squeeze(sum(var_signal_k,1));
        var_sum = sum(var_signal_k(:));
        az_spread(row,col) = P.var_theta/P.dtheta^2*az_var_signal/var_sum;
        rad_spread(row,col) = P.var_rad/P.drad^2*rad_var_signal'/var_sum;
        
        k = k + 1;
    end
end

%% Plot 
figure(1)
imshow(az_spread,'DisplayRange',[1 256],'Colormap',jet,'InitialMagnification','fit')
figure(2)
imshow(rad_spread,'DisplayRange',[1 4],'Colormap',jet,'InitialMagnification','fit')
figure(3)
imshow(rel_err,'DisplayRange',[0 1],'Colormap',jet,'InitialMagnification','fit')
figure(4)
subplot(2,1,1)
imshow(img_fit,'DisplayRange',[0 max(img_fit(:))],'Colormap',jet,'InitialMagnification','fit')
subplot(2,1,2)
imshow(polar_image/norm(polar_image(:)),'DisplayRange',[0 max(img_fit(:))],'Colormap',jet,'InitialMagnification','fit')