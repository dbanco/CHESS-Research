function out = ToyTest5x5()
%ToyTest5x5 Summary of this function goes here
%   Detailed explanation goes here

T = 29;
s = zeros(9,9,T);
i = 3;
j = 3;
steps = [1 1; -1 0; 0 -1; 1 0; 0 1;-1 -1; 0 0; 0 0];
set = 1;
setSize = 4;

for t = 1:T
    s(i,j,t) = 1;
%     s(i,j+1,t) = 1;
%     s(i+1,j+1,t) = 1;
    s(i+1,j,t) = 1;
    s(i-1,j,t) = 1;
    i = i + steps(set,1);   
    j = j + steps(set,2);
    if ~mod(t,setSize)
        set = set + 1;
    end
end

blobS = convn(s,ones(3,3),'same');

smoothness = 1;
maxIters = 10;
[u,v,Fy,Fx,Ft] = computeHornSchunk(s,smoothness,maxIters);
opticFlow = opticalFlowHS('Smoothness',smoothness,'MaxIteration',maxIters);
term = s;
term(Ft == 0) = 0;
of_con = opticalFlowOp(s,u,v,1,0,0);
out = norm(vec(of_con))
of_con2 = opticalFlowOp(s,u,v,1,0,1);
out2 = norm(vec(of_con2))


%% View dataset
close all

h1 = figure;
h1.Position = [2.5826e+03 391.4000 560 420.0000];

h2 = figure;
h2.Position = [2.5942e+03 -108.6000 560 420]; 

h11 = figure;
h11.Position = [3.1534e+03 387.8000 560 420.0000];
h12 = figure;
h12.Position = [3.7294e+03 381 560.0000 420];

h21 = figure;
h21.Position = [3.1662e+03 -113 560 420];
h22 = figure;
h22.Position = [3.7322e+03 -118.2000 560.0000 420];

h3 = figure;
h3.Position = [4375 386.6000 560 420];

h4 = figure;
h4.Position = [4375 386.6000 560 420];

h5 = figure;
h5.Position = [4375 386.6000 560 420];

h6 = figure;
h6.Position = [4375 386.6000 560 420];




for t = 1:T
    frameGray = squeeze(s(:,:,t));
    figure(h1)
    flow = estimateFlow(opticFlow,frameGray);
    imagesc(frameGray)
    hold on
    quiver(flow.Vx,flow.Vy)
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off
    title('matlab')

    figure(h2)
    imagesc(frameGray)
    hold on
    quiver(u(:,:,t),v(:,:,t))
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off
    title(sprintf('I(\\eta,\\sigma,%i)',t))

    figure(h11)
    imagesc(flow.Vx)
    title('Vx')

    figure(h12)
    imagesc(flow.Vy)
    title('Vy')

    figure(h21)
    imagesc(u(:,:,t))
    title('u')

    figure(h22)
    imagesc(v(:,:,t))
    title('v')
    

    figure(h3)
    imagesc(of_con(:,:,t))
    title('optical flow constraint')

%     figure(h4)
%     imagesc(Ft(:,:,t))
%     title('F_t')
% 
%     figure(h5)
%     imagesc(Fx(:,:,t))
%     title('F_x')
% 
%     figure(h6)
%     imagesc(Fy(:,:,t))
%     title('F_y')
    pause()
end
%% Write gif
gifDir = 'C:\Users\dpqb1\Desktop\HStoy_ex\c';
fName = 'flowHSc_1_10.gif';
h2 = figure;
h2.Position = [2.5942e+03 -108.6000 560 420];
for t = 1:T
    figure(h2)
    frameGray = squeeze(s(:,:,t));
    imagesc(frameGray)
    hold on
    quiver(u(:,:,t),v(:,:,t))
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off
    title(sprintf('I(\\eta,\\sigma,%i)',t))
    
    pause(0.001)

    frame = getframe(h2);
    im = frame2im(frame);

    [imind,cm] = rgb2ind(im,256);

    % Save the image to a GIF file
    if t == 1
        outFile = fullfile(gifDir,fName);
        imwrite(imind,cm,outFile,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,outFile,'gif','WriteMode','append','DelayTime',0.2);
    end
end

%% 
% t = 2;
% 
% figure
% subplot(5,1,1)
% imagesc(Fx(:,:,t))
% title('F_x')
% subplot(5,1,2)
% imagesc(u(:,:,t))
% title('u')
% subplot(5,1,3)
% imagesc(Fy(:,:,t))
% title('F_y')
% subplot(5,1,4)
% imagesc(v(:,:,t))
% title('v')
% subplot(5,1,5)
% imagesc(Ft(:,:,t))
% title('F_t')
% 
% 
% figure
% subplot(5,1,1)
% imagesc(Fx(:,:,t).*u(:,:,t))
% title('u*F_x')
% subplot(5,1,2)
% imagesc(Fy(:,:,t).*v(:,:,t))
% title('v*F_y')
% subplot(5,1,3)
% imagesc(Ft(:,:,t)+Fy(:,:,t).*v(:,:,t)+Fx(:,:,t).*u(:,:,t))
% title('OF')
% subplot(5,1,4)
% imagesc(s(:,:,t))
% title('s_{t}')
% subplot(5,1,5)
% imagesc(s(:,:,t))
% title('s_{t+1}')


