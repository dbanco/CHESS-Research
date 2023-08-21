function plotOpticalFlow(X,K,opt,fName,outputDir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin > 3
    gifDir = outputDir;
    outFile = fullfile(gifDir,fName);
end
X = squeeze(X);
[N,KJ,T] = size(X);
J = KJ/K;
% Compute opticalflow
[u,v] = computeHornSchunkDict(X,K,opt.Smoothness,opt.HSiters);

% Show velocity in time
% avgU = mean(mean(u,1),2);
% avgV = mean(mean(v,1),2);
% timeU = zeros(N,KJ);
% timeV = zeros(N,KJ);
% % Find max of velocity
% figure
% for t = 1:100
%     [valU,ui] = max(u(:,:,t),[],'all');
%     [valV,vi] = max(v(:,:,t),[],'all');
%     timeU(ui) = timeU(ui) + valU;
%     timeV(ui) = timeV(ui) + valV;
% 
% end
% quiver(timeV',timeU')
% hold on
% Compute optical flow constraint
% OF = sum( (opticalFlowOp(X,u,v,K)).^2,'all')

window = 50:100;

% View optical flow
fig = figure;
fig.Position = [1 1 1.6452e+03 554.8000];
opticFlow1 = opticalFlowHS;
opticFlow2 = opticalFlowFarneback;
opticFlow3 = opticalFlowLK;
for t = 1:T
    frameGray = squeeze(X(:,1:J,t));
    flow1 = estimateFlow(opticFlow1,frameGray);
    flow2 = estimateFlow(opticFlow2,frameGray);
    flow3 = estimateFlow(opticFlow3,frameGray);

    subplot(4,1,1)
    framePlot = frameGray(window,:)';
    imagesc(framePlot)
    hold on
    quiver(v(window,1:J,t)',u(window,1:J,t)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off

    subplot(4,1,2)
    imagesc(framePlot)
    hold on
    quiver(flow1.Vy(window,:)',flow1.Vx(window,:)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off

    subplot(4,1,3)
    imagesc(framePlot)
    hold on
    quiver(flow2.Vy(window,:)',flow2.Vx(window,:)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off

    subplot(4,1,4)
    imagesc(framePlot)
    hold on
    quiver(flow3.Vy(window,:)',flow3.Vx(window,:)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off
    % Save the image to a GIF file
    if nargin > 3
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if t == 1
            imwrite(imind,cm,outFile,'gif','Loopcount',inf,'DelayTime',0.2);
        else
            imwrite(imind,cm,outFile,'gif','WriteMode','append','DelayTime',0.2);
        end
    end
    pause()
end

end