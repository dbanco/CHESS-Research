function plotOpticalFlow4param(X,K,opt,fName,outputDir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin > 3
    gifDir = outputDir;
    outFile = fullfile(gifDir,fName);
end
X = squeeze(X);
[N,KJ,T] = size(X);
J = KJ/K;
window = 1:size(X,1);

% Compute opticalflow
HS_params = [1e-8,1e-4,1e-2,1,10];
L = numel(HS_params);
u_array = cell(L,1);
v_array = cell(L,1);
for i = 1:L
    [u,v] = computeHornSchunkDict(X,K,HS_params(i),opt.HSiters);
    u_array{i} = u;
    v_array{i} = v;
end

% View optical flow
fig = figure;
fig.Position = [1 1 1400 900];

for t = 1:4
    frameGray = squeeze(X(:,:,t));

    framePlot = frameGray(window,:)';
    
    for i = 1:L
        u = u_array{i};
        v = v_array{i};
        
        u = u.*X;
        v = v.*X;

        u = u/max(u(:))*100;
        v = v/max(v(:))*100;

        subplot(L,1,i)
        imagesc(framePlot)
        hold on
        quiver(v(window,:,t)',u(window,:,t)')
        q = findobj(gca,'type','Quiver');
        q.Color = 'r';
        hold off
    end

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