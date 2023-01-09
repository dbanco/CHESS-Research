% view omega behavior of mmpad data
data_dir = 'E:\\MMPAD_omega\\full_etaSum\\';
% Notes: Weird artifact at omega = 34, t = 47
r1 = [5,210,268,355]; r2 = [40,260,310,395];
title_str = {'\{004\}','\{021\}','\{112\}','\{020\}'};
sweeps = [1,73
          8,65;...
          15,58;...
          22,50];

% Summing in eta to look only at omega in time
fname = ['mmpad_img_',num2str(1),'.mat'];
load(fullfile(data_dir,fname))
[N1,N2,N3] = size(mmpad_img);
T = 300;
R = 4;
B = zeros(N1,T,R);
for i = 1:T
    if ~mod(i,10)
        fprintf('%i\n',i)
    end
    fname = ['mmpad_img_',num2str(i),'.mat'];
    load(fullfile(data_dir,fname))
    for r = 1:R
        B(:,i,r) = sum(mmpad_img(:,r1(r):r2(r)),2);
    end
end
f1 = figure(1);
for r = 4
    imagesc(B(:,:,r)')
    xlabel('\omega')
    ylabel('time')
    title(title_str{r})
    hold on
    for i = 1:4
        plot([sweeps(i,1),sweeps(i,1)],[1,T],'r')
        plot([sweeps(i,2),sweeps(i,2)],[1,T],'r')
    end
end

f1.Position = [971,573,560,420];

%% Inspect data
%{
for i = 0:200
    fname = [' mmpad_img_',num2str(i),'.mat'];
    gifname = 'eta_omega_mmpad_021';
    load(fullfile(data_dir,fname))

    f1 = figure(1);
    imagesc(squeeze(sum(mmpad_img(:,:,r1:r2),3)))
    ylabel('\omega')
    xlabel('\eta')
    title(['t = ',num2str(i)])
    f1.Position = [971,573,560,420];
    
    % creat gif
    drawnow
    frame = getframe(f1);
    [A,map] = rgb2ind(frame2im(frame),256);
    if i==0
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.1);
    end
%     f2 = figure(2);
%     imagesc(squeeze(sum(mmpad_img(:,:,r1:r2),1))')
%     ylabel('2\theta')
%     xlabel('\eta')  
%     title(['t = ',num2str(i)])
%     f2.Position = [968,408,560,81];
end
%}