% view omega behavior of mmpad data
data_dir = 'E:\\MMPAD_omega\\full\\';
% Notes: Weird artifact at omega = 34, t = 47

%% 004 Ring
r1 = 5; r2 = 40;
title_str = '\{004\}';
subDir = 'ring1';
% Notes:
% Spot exits at eta = 80, omega = 1, t = 20
% Spot enters at omega = 48, eta = 1, t= 21
% Some signal enters at eta = 238, omega = 1, t = 70
% (artifact t=47 visible)
%% 021 Ring
r1 = 210; r2 = 260;
title_str = '\{021\}';
subDir = 'ring1';
% Notes:
% Spot enters omega = 53, eta = 1, t = 25
% Spot enters eta = 193, omega = 72, t = 29
% Spot partial omega = 40, eta = 265, t = 42
% Big signal exits eta = 139, t = 79
% Big signal exits eta = 1, t = 79
% Signal entering eta = 72, t = 95
% (artifact t=47 visible)
%% 112 Ring
r1 = 268; r2 = 310;
title_str = '\{112\}';
% Notes:
% Spot enters eta = 17, omega = 72, t = 21
% Spots exit eta = 100, omega = 1, t = 32
% Spot exits omega = 22, eta = 265, t = 38
% Big signal exits eta = 19, omega = 72, t = 66
% Signal enters eta = 171, omega = 1, t = 112
% (artifact t=47 not visible)
%% 020 Ring
r1 = 355; r2 = 395;
title_str = '\{020\}';
% Notes:
% Spot exits eta = 219, omega = 1, t = 23
% Spot enters eta = 55, omega = 72, t = 32
% Spot exits eta = 251, omega = 1, t = 33
% (artifact t=47 visible)

%% Summing in eta to look only at omega in time
fname = ['mmpad_img_',num2str(1),'.mat'];
load(fullfile(data_dir,fname))
[N1,N2,N3] = size(mmpad_img);
T = 546;
B = zeros(N1,T);
for i = 0:T-1
    if ~mod(i,10)
        fprintf('%i\n',i)
    end
    fname = ['mmpad_img_',num2str(i),'.mat'];
    load(fullfile(data_dir,fname))
    mmpad_img = squeeze(sum(mmpad_img,2));
    save(fullfile('E:\MMPAD_omega\full_etaSum',...
                  ['mmpad_img_',num2str(i+1),'.mat']),'mmpad_img')
    B(:,i+1) = sum(mmpad_img(:,r1:r2),2);
end
f1 = figure(1);
imagesc(B)
ylabel('\omega')
xlabel('time')
title(title_str)
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