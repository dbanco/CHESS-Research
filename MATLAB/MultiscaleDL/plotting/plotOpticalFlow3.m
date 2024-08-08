function plotOpticalFlow3(X,K,opt,fName,outputDir)
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
[u1,v1,~,~,~] = computeHornSchunkDict(X,K,opt.Smoothness,opt.HSiters);
[u2,v2,~,~,~] = computeHornSchunkDictLS(X,K,[],[],opt.Smoothness,opt.HSiters);
[u3,v3,~,~,~] = computeHornSchunkDictPaper(X,K,opt.Smoothness,opt.HSiters);
[u4,v4,~,~,~] = computeHornSchunkDictPaperLS(X,K,opt.Smoothness,opt.HSiters);

window =1:300;

% View optical flow
fig = figure;
fig.Position = [1 1 1.6452e+03 554.8000];
for t = 1:T
    frameGray = squeeze(X(:,1:J,t));
    framePlot = frameGray(window,:)';

    subplot(4,1,1)
    imagesc(framePlot)
    hold on
    quiver(v1(window,1:J,t)',u1(window,1:J,t)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off

    subplot(4,1,2)
    imagesc(framePlot)
    hold on
    quiver(v2(window,1:J,t)',u2(window,1:J,t)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off

    subplot(4,1,3)
    imagesc(framePlot)
    hold on
    quiver(v3(window,1:J,t)',u3(window,1:J,t)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off

    subplot(4,1,3)
    imagesc(framePlot)
    hold on
    quiver(v4(window,1:J,t)',u4(window,1:J,t)')
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