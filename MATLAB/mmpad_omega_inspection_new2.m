% view omega behavior of mmpad data
data_dir = 'E:\\MMPAD_omega';
% Notes: Weird artifact at omega = 34, t = 47
r1 = [5,210,268,355]; r2 = [40,260,310,395];
title_str = {'\{004\}','\{021\}','\{112\}','\{020\}'};
omDir = {'omega2','omega3','omega4','omega5'};
rDir = {'ring1','ring2','ring3','ring4'};
sweeps = [1,73
          8,65;...
          15,58;...
          22,50];

% Summing in eta to look only at omega in time
fname = ['mmpad_img_',num2str(1),'.mat'];
load(fullfile(data_dir,omDir{1},rDir{1},fname))
[N1,~] = size(polar_image);
T = 300;
R = 4;
Oms = 4;
B = zeros(N1,T,Oms,R);
for r = 1:4
    fprintf('Ring %i Omega ',r)
    for o = 1:Oms
        fprintf('%i, ',o)
        for t = 1:T
            if ~mod(t,10)
                
            end
            fname = ['mmpad_img_',num2str(t),'.mat'];
            load(fullfile(data_dir,omDir{o},rDir{r},fname))
           B(:,t,o,r) = sum(polar_image,2);
        end
    end
    fprintf('\n')
end

%%
f1 = figure(1);
for r = 3
    for o = 1:4
        subplot(1,4,o)
        imagesc( (B(:,1:80,o,r)./max(B(:,1:80,o,r)))' )
        xlabel('\eta')
        ylabel('time')
        title(omDir{o})
    end
end

f1.Position = [971,573,560,420];

figure(2)
plot(B(:,50,4,3))

figure(3)
plot(B(:,37,4,3))

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