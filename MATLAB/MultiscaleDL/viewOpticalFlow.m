% Load outputs
outputDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF4_sig_2';
load(fullfile(outputDir,'output_1.mat'));
data = squeeze(outputs.Xmin);
[N, KJ, T] = size(data);
K = 2;

% Compute opticalflow
[u,v] = computeHornSchunkDict(data,2);

% Compute optical flow constraint
OF = sum( (opticalFlowOp(data,u,v,K)).^2,'all')

% View optical flow
fig = figure;
fig.Position = [2716.2 313.8 1.6452e+03 554.8000];
gifDir = outputDir;
fName = 'opticalFlow_1min.gif';
opticFlow = opticalFlowHS;
for t = 1:T
    frameGray = squeeze(data(:,:,t));
    imagesc(frameGray')
    hold on
    quiver(v(:,:,t)',u(:,:,t)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off
    
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Save the image to a GIF file
    if t == 1
        outFile = fullfile(gifDir,fName);
        imwrite(imind,cm,outFile,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,outFile,'gif','WriteMode','append','DelayTime',0.2);
    end

    pause()
end