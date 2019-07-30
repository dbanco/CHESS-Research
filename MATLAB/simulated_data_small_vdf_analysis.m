%% Show video
video_norm = cell(num_imgs,1);

baseDir = 'D:\CHESS_data\';
datasetName = 'two_spot_growth_init';
fitName = '';
fDir = [baseDir,datasetName,fitName];

for img_num = 1:num_imgs
    frame = [];
    fprintf('Image %i\n',img_num)
    fName = sprintf('fista_fit_%i_%i.mat',lam_num,img_num);
    load(fullfile(fDir,fName))
    
    A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
    img_fit = Ax_ft_2D(A0ft_stack,x_hat)*norm(polar_image(:));
    err_fit = norm(img_fit-polar_image)/norm(polar_image);
    vdf = squeeze(sum(sum(x_hat,1),2))/sum(x_hat(:));
    % Append polar_images and fits
    frame = [frame,polar_image];
    video_norm{img_num} = frame;
        lim1 = 0;
        lim2 = max(polar_image(:));
        lim_vdf = max(vdf(:));
        % Plot both images
        figure(11)
        subplot(1,3,1)
        imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
        subplot(1,3,2)
        imshow(img_fit,'DisplayRange',[lim1 lim2],'Colormap',jet);
        title(sprintf('Error = %f',err_fit))
        subplot(1,3,3)
        imshow(vdf,'DisplayRange',[0 lim_vdf],'Colormap',jet);
        title(VDF)
        pause(0.1)
end

%% Create gif
h = figure;
filename = [datasetName,fitName,'.gif'];
for i = 1:num_imgs
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

    [imind,cm] = rgb2ind(im_crop,256,'nodither'); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
    end 
end

%% Show mmpad vdf video
video_norm_diff = cell(100,1);

fDir = ['D:\CHESS_data\',datasetName,fitName];
for img_num = 1:100
    display(img_num)
    frame = [];
    for ring_num = 1
        fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);
        
        load(fullfile(fDir,fName))
%         A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
%         img_fit = Ax_ft_2D(A0ft_stack,x_hat)*norm(polar_image(:));
        vdf1 = squeeze(sum(sum(x_hat,1),2));
        % Append polar_images and fits
        frame = abs(vdf1);
        vdf_i = VDF{img_num};
        vdf_i = vdf_i/sum(vdf_i(:));
        frameTrue = abs(vdf_i);
    end
    video_norm_diff{img_num} = frame;
    video_norm_true{img_num} = frameTrue;
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
filename = [datasetName,fitName,'_vdf.gif'];
for i = 1:100
    lim1 = 0;
    lim2 = max(video_norm_diff{i}(:));
    imshow(video_norm_diff{i},'DisplayRange',[lim1 lim2],'Colormap',jet);
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
    pixels = 2;
    
    ring_num = 1;
    
    [imind,cm] = rgb2ind(im_crop,256,'nodither'); 
    % Write to the GIF File 
    if i == 1
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
    end 
end