baseDir = 'D:\CHESS_data\wass_small_fit\';
datasetName = 'wass_small';
fitName = '_fit_b';

num_imgs = 20;
lambda_vals = [0.00001 0.00002 0.00005 0.0001 0.0002 0.0005  0.001  0.002,...
                0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5];
% Load fit
img_num = 10;
ring_num = 1;

fDir = [baseDir,datasetName,fitName];
fName = sprintf('fista_fit_%i_%i.mat',ring_num,img_num);

load(fullfile([fDir,'_1'],fName))
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


% Construct distance matrix
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


% Compute time varying spread and error functions
P.num_var_t = 15;
P.num_var_r = 10;

fDir = [baseDir,datasetName,fitName];
az_spread = zeros(num_imgs,numel(lambda_vals));
rad_spread = zeros(num_imgs,numel(lambda_vals));
rel_err = zeros(num_imgs,numel(lambda_vals));
sparsity = zeros(num_imgs,numel(lambda_vals));
%wass_dist = zeros(num_imgs,1);
var_signal = zeros(num_imgs,numel(lambda_vals),P.num_var_t,P.num_var_r);

for lam_num = 1:numel(lambda_vals)
    for img_num = 1:num_imgs
        k = img_num;
        fprintf('Image %i\n',k)
        fName = sprintf('fista_fit_%i_%i.mat',lam_num,img_num);
        load(fullfile([fDir,'_',num2str(lam_num)],fName))
        A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
        img_fit = Ax_ft_2D(A0ft_stack,x_hat);
        rel_err(k,lam_num) = err(end);
        var_signal_k = squeeze(sum(sum(x_hat,1),2));
        var_signal(k,lam_num,:,:) =  var_signal_k;
        rad_var_signal = squeeze(sum(var_signal_k,2));
        az_var_signal = squeeze(sum(var_signal_k,1));
        var_sum = sum(var_signal_k(:));
        rad_spread(k,lam_num) = sqrt(P.var_rad)*az_var_signal'/var_sum;
        az_spread(k,lam_num) = sqrt(P.var_theta)*rad_var_signal/var_sum;
        sparsity(k,lam_num) = sum(x_hat(:)>0);

    %     if k > 2
    %         wass_dist(k) = sinkhornKnoppTransport(var_signal_k(:),vdf_last(:),P.params.wLam,D);
    %     end

        k = k + 1;
        vdf_last = var_signal_k;
    end
end
    spreadDir = fullfile('D:','CHESS_data','spread_results');
    mkdir(spreadDir)
    outFile = ['spread_',datasetName,fitName,'.mat'];
    save(fullfile(spreadDir,sprintf(outFile,ring_num)),...
        'var_signal','rel_err','P','rad_spread','az_spread','rel_err','sparsity')%,'wass_dist')
 
%% Load spread data
spreadDir = fullfile('D:','CHESS_data','spread_results');

ring_data{1} = load(fullfile(spreadDir,['spread_',datasetName,fitName,'.mat']));

load([baseDir,'wass_small','\synth_data.mat'])
truth_awmv_az = zeros(num_imgs,1);
truth_awmv_rad = zeros(num_imgs,1);
for i = 1:num_imgs
    vdf_i = VDF{i};
    vdf_i = vdf_i/sum(vdf_i(:));
    truth_awmv_az(i) = sum(vdf_i,2)'*sqrt(P.var_theta)';
    truth_awmv_rad(i) = sum(vdf_i,1)*sqrt(P.var_rad)';
end
%% Plot time varying spread and error functions
close all
colors = jet(numel(lambda_vals)+1);

num_lines = 18;

figure(1)
for i = 1:num_lines
    cv = colors(i,:);
    plot(ring_data{1}.rad_spread(:,i),'-o','Color',cv,'MarkerSize',3)
    hold on
end
plot(truth_awmv_rad,'-x','MarkerSize',2)
leg_str = cell(1,numel(lambda_vals));
for i = 1:numel(lambda_vals)
   leg_str{1,i} = num2str(lambda_vals(i));
end
leg_str{1,end+1} = 'truth';
legend(leg_str,'Location','Best')
title('Radial AWMV')
xlabel('Time')

legend_vals = ['0.00001', '0.00002', '0.00005', '0.0001','0.0002', '0.0005',  '0.001', '0.002',...
                '0.005', '0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1', '2', '5'];

figure(2)
for i = 1:num_lines
    cv = colors(i,:);
    plot(ring_data{1}.az_spread(:,i),'-o','Color',cv,'MarkerSize',3)
    hold on
end
plot(truth_awmv_az,'-x','MarkerSize',2)
leg_str = cell(1,numel(lambda_vals));
for i = 1:numel(lambda_vals)
   leg_str{1,i} = num2str(lambda_vals(i));
end
leg_str{1,end+1} = 'truth';
legend(leg_str,'Location','Best')
title('Azimuthal AWMV')
xlabel('Time')

figure(3)
for i = 1:num_lines
    cv = colors(i,:);
    plot(ring_data{1}.rel_err(:,i),'-o','Color',cv,'MarkerSize',3)
    hold on
end
leg_str = cell(1,numel(lambda_vals));
for i = 1:numel(lambda_vals)
   leg_str{1,i} = num2str(lambda_vals(i));
end
legend(leg_str,'Location','Best')
title('Relative Error')
xlabel('Time')

figure(4)
for i = 1:num_lines
    cv = colors(i,:);
    semilogy(ring_data{1}.sparsity(:,i),'-o','Color',cv)
    hold on
end
leg_str = cell(1,numel(lambda_vals));
for i = 1:numel(lambda_vals)
   leg_str{1,i} = num2str(lambda_vals(i));
end
legend(leg_str,'Location','Best')
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
video_norm = cell(num_imgs,1);
fDir = ['D:\CHESS_data\',datasetName,fitName,'_1'];

for img_num = 1:num_imgs
    display(img_num)
    frame = [];
    fName = sprintf('fista_fit_%i_%i.mat',1,img_num);
    load(fullfile(fDir,fName))
    A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
    img_fit = Ax_ft_2D(A0ft_stack,x_hat)*norm(polar_image(:));

    % Append polar_images and fits
    frame = [frame,polar_image];
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