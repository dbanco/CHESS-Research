%% Show mmpad video
video_norm = cell(546,1);

for img_num = 1:546
    display(img_num)
    frame = [];
    for ring_num = 1:4
        fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);
        fDir = ['D:\MMPAD_data\mmpad_ring',num2str(ring_num),'_zero_fit_reg4'];
        load(fullfile(fDir,fName))
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
        img_fit = Ax_ft_2D(A0ft_stack,x_hat)*norm(polar_image(:));
        
        % Append polar_images and fits
        frame = [frame,polar_image,img_fit];
    end
    video_norm{img_num} = frame;
%         lim1 = 0;
%         lim2 = max(polar_image(:));
%         % Plot both images
%         figure(11)
%         subplot(1,2,1)
%         imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
%         subplot(1,2,2)
%         imshow(img_fit,'DisplayRange',[lim1 lim2],'Colormap',jet);
%         title(sprintf('Error = %f',err(end)))
%         pause
end

%% Create gif
h = figure;
filename = 'mmpad_norm_fit4.gif';
for i = 1:546
    lim1 = 0;
    lim2 = max(video_norm{i}(:));
    imshow(video_norm{i},'DisplayRange',[lim1 lim2],'Colormap',jet);
    drawnow
    frameg = getframe(h);
    im = frame2im(frameg);
    % Remove grayspace to crop image
    m_im = mean(im,3);
    im_rows = mean(m_im,1);
    im_cols = mean(m_im,2);
    rows = find(im_cols~=240);
    cols = find(im_rows~=240);
    im_crop = im(rows,cols,:);
    [nn,mm,pp] = size(im_crop);
    % Add whitespace between rings and whiten non detected pixels
    ind = 0;
    pixels = [2,10];
    for ring_num = 1:4
        for dd = 1:2
            fName = sprintf('fista_fit_%i_%i.mat',ring_num,1);
            fDir = ['D:\MMPAD_data\mmpad_ring',num2str(ring_num),'_zero_fit_reg4'];
            load(fullfile(fDir,fName))
            n_cols = size(polar_image,2);
            ind = ind + n_cols;
            if(dd==1)
                im_crop(129:133,ind-n_cols+1:ind,:) = 256;
            end
            im_crop = [im_crop(:,1:ind,:),256*ones(nn,pixels(dd),pp),im_crop(:,ind+1:end,:)];
            ind = ind + pixels(dd);
        end
    end
    [imind,cm] = rgb2ind(im_crop,256); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end

%% Load mmpad fit
img_num = 25;
ring_num = 1;

fDir = 'D:\MMPAD_data\mmpad_ring1_zero_fit_reg5_init';
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
imshow(img_fit*norm(polar_image(:)),'DisplayRange',[lim1 lim2],'Colormap',jet);

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
for ring_num = 1
    fDir = ['D:\MMPAD_data\mmpad_ring',sprintf('%i',ring_num),'_zero_fit_reg5_init'];
    az_spread = zeros(546,1);
    rad_spread = zeros(546,1);
    rel_err = zeros(546,1);
    var_signal = zeros(546,P.num_var_t,P.num_var_r);
    k = 1;
    for img_num = 1:546
        fprintf('Image %i\n',k)
        fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);
        load(fullfile(fDir,fName))
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

    spreadDir = fullfile('D:','MMPAD_data','spread_results');
    mkdir(spreadDir)
    outFile = 'spread_mmpad_ring%i_zero_fit5_init_lambda_%2.3f.mat';
    save(fullfile(spreadDir,sprintf(outFile,ring_num,P.params.lambda)),...
        'var_signal','rel_err','P','rad_spread','az_spread','rel_err')
end    
%% Load spread data
for i = 1
    ring_data{i} = load(fullfile(spreadDir,sprintf('spread_mmpad_ring%i_zero_fit5_init_lambda_0.010.mat',i)));
end
%% Plot time varying spread and error functions

figure(1)
for j = 1
    plot(ring_data{j}.rad_spread(1:200),'-o','MarkerSize',2)
    hold on
end
legend('1','2','3','4','Location','Best')
title('Azimuthal AWMV')
xlabel('Time')

figure(2)
for j = 1
    plot(ring_data{j}.az_spread(1:200),'-o','MarkerSize',2)
    hold on
end
legend('1','2','3','4','Location','Best')
title('Radial AWMV')
xlabel('Time')


figure(3)
for j = 1
    plot(ring_data{j}.rel_err,'-')
    hold on
end
legend('1','2','3','4','Location','Best')
title('Relative Error')
xlabel('Time')

%% Show evolving var_signal
lim1 = 0;
lim2 = max(var_signal(:));
for i = 1:546
    imshow(squeeze(var_signal(i,:,:)),'DisplayRange',[lim1 lim2],'Colormap',jet,'InitialMagnification','fit');
    pause(0.2)
end