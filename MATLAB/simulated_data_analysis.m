datasetName = 'simulated_data2';
fitName = '_fit_wass_reg';

%% Load fit
img_num = 10;
ring_num = 1;

fDir = ['D:\CHESS_data\',datasetName,fitName];
fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);

load(fullfile(fDir,fName))
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
img_fit = Ax_ft_2D(A0ft_stack,x_hat);

lim1 = 0;
lim2 = max(polar_image(:));
% Plot both images
% figure(1)
% subplot(2,1,1)
% imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
% subplot(2,1,2)
% imshow(img_fit*norm(polar_image(:)),'DisplayRange',[lim1 lim2],'Colormap',jet);

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

%% Construct distance matrix
N = P.num_var_t*P.num_var_r;
THRESHOLD = 32;

D = ones(N,N).*THRESHOLD;
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        for ii=max([1 i-THRESHOLD+1]):min([P.num_var_t i+THRESHOLD-1])
            for jj = max([1 j-THRESHOLD+1]):min([P.num_var_r j+THRESHOLD-1])
                ind1 = i + (j-1)*P.num_var_t;
                ind2 = ii + (jj-1)*P.num_var_t;
                D(ind1,ind2)= sqrt((i-ii)^2+(j-jj)^2); 
            end
        end
    end
end
D = D./max(D(:));


%% Compute time varying spread and error functions
P.num_var_t = 15;
P.num_var_r = 10;

fDir = ['D:\CHESS_data\',datasetName,fitName];
az_spread = zeros(100,1);
rad_spread = zeros(100,1);
rel_err = zeros(100,1);
sparsity = zeros(100,1);
%wass_dist = zeros(100,1);
var_signal = zeros(100,P.num_var_t,P.num_var_r);

for img_num = 1:100
    k = img_num;
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
    sparsity(k) = sum(x_hat(:)>0);

%     if k > 2
%         wass_dist(k) = sinkhornKnoppTransport(var_signal_k(:),vdf_last(:),P.params.wLam,D);
%     end

    k = k + 1;
    vdf_last = var_signal_k;
end
    spreadDir = fullfile('D:','CHESS_data','spread_results');
    mkdir(spreadDir)
    outFile = ['spread_',datasetName,fitName,'.mat'];
    save(fullfile(spreadDir,sprintf(outFile,ring_num)),...
        'var_signal','rel_err','P','rad_spread','az_spread','rel_err','sparsity')%,'wass_dist')
 
%% Load spread data
spreadDir = fullfile('D:','CHESS_data','spread_results');

ring_data{1} = load(fullfile(spreadDir,['spread_',datasetName,fitName,'.mat']));

load(['D:\CHESS_data\','simulated_data1','\synth_data.mat'])
truth_awmv_az = zeros(100,1);
truth_awmv_rad = zeros(100,1);
for i = 1:100
    vdf_i = VDF{i};
    vdf_i = vdf_i/sum(vdf_i(:));
    truth_awmv_rad(i) = sum(vdf_i,2)'*sqrt(P.var_theta)';
    truth_awmv_az(i) = sum(vdf_i,1)*sqrt(P.var_rad)';
end
%% Plot time varying spread and error functions
figure(1)
for j = 1
    plot(ring_data{j}.rad_spread(1:100),'-o','MarkerSize',3)
    hold on
    plot(truth_awmv_rad,'-x','MarkerSize',2)
end
legend('Fit','Truth','3','4','Location','Best')
title('Azimuthal AWMV')
xlabel('Time')

figure(2)
for j = 1
    plot(ring_data{j}.az_spread(1:100),'-o','MarkerSize',3)
    hold on
    plot(truth_awmv_az,'-x','MarkerSize',2)
end
legend('Fit','Truth','3','4','Location','Best')
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

figure(4)
for j = 1
    semilogy(ring_data{j}.sparsity,'o')
    hold on
end
%legend(num2str(mean(ring_data{1}.sparsity)),num2str(mean(ring_data{2}.sparsity)),'Location','Best')
title('Sparsity')
xlabel('Time')

% figure(6)
% for j = 1
%     plot(ring_data{j}.wass_dist,'-')
%     hold on
% end
% legend('1','2','3','4','Location','Best')
% title('Wasserstein distance to neighbor')
% xlabel('Time')

%% Show evolving var_signal
figure()
% Load truth vdf
% load(['D:\CHESS_data\',datasetName,'\synth_data.mat'])
for i = 1:100
    im_vdf = video_norm_diff{i};
    im_vdf = im_vdf/sum(im_vdf(:));
    lim1 = min(im_vdf(:));
    lim2 = max(im_vdf(:));
    subplot(1,3,1)
    imshow(im_vdf,'DisplayRange',[lim1 lim2],'Colormap',jet,'InitialMagnification','fit');
    
    % Compute matching vdf
    [ind_theta,ind_rad] = find(ones(size(im_vdf)));
    mean_theta = sum(im_vdf(:).*ind_theta);
    mean_rad = sum(im_vdf(:).*ind_rad);
    vdf_match = gaussian_basis_2D(size(im_vdf,1),mean_theta,0.25, size(im_vdf,2),mean_rad,0.25);
    vdf_match = reshape(vdf_match,size(im_vdf));
    subplot(1,3,2)
    imshow(vdf_match,'DisplayRange',[lim1 lim2],'Colormap',jet,'InitialMagnification','fit');
    
    % Display true vdf
    vdf_true = VDF{i};
    vdf_true = vdf_true/sum(vdf_true(:));
    subplot(1,3,3)
    imshow(vdf_true,'DisplayRange',[lim1 lim2],'Colormap',jet,'InitialMagnification','fit');
    
    pause(0.2)
end


%% Show mmpad video
video_norm = cell(100,1);
fDir = ['D:\CHESS_data\',datasetName,fitName];

for img_num = 1:100
    display(img_num)
    frame = [];
    fName = sprintf('fista_fit_%i_%i.mat',1,img_num);
    load(fullfile(fDir,fName))
    A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
    img_fit = Ax_ft_2D(A0ft_stack,x_hat)*norm(polar_image(:));

    % Append polar_images and fits
    frame = [frame,polar_image,img_fit];
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
filename = [datasetName,fitName,'.gif'];
for i = 1:100
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