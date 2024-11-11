function plotDataRecon(y,recon,gifDir,fName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y = squeeze(y);
recon = squeeze(recon);
T = size(y,2);
yMax = max(y(:));

fig = figure;
for t = 1:T
    plot(y(:,t))
    hold on
    plot(recon(:,t))
    hold off
    
    ylim([0 yMax*1.1])
%         pause()
    if nargin > 1
        % Capture the current frame as an image
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        % Save the image to a GIF file
        if t == 1
            outFile = fullfile(gifDir,fName);
            imwrite(imind,cm,outFile,'gif','Loopcount',inf,'DelayTime',0.1);
        else
            imwrite(imind,cm,outFile,'gif','WriteMode','append','DelayTime',0.1);
        end
    else
%         pause()
    end
end



end

