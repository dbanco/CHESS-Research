% Optical flow example
data = gaus_linear_osc_signal_2D(0.002);
% data = loadMMPAD2D(1,1,'C:\Users\dpqb1\Documents\Data\ring1_zero');
% load('C:\Users\dpqb1\Documents\Outputs\toy1_exp1_sig_1\output_22.mat')
% data = squeeze(outputs.X);
[Nx,Ny,Nt] = size(data);

%% Test Transpose operators
tic
Fy = sobel(data,1);
toc
tic
Fty = sobel(data,1,1);
toc

% Fx = sobel(data,2);
% Ft = -timeDiff(data);
% 
% Fty = sobel(data,1);
% Ftx = sobel(data,2);

%% Use Lucas and Kanade method to compute optical flow locally
% win = 6;
% V = localFlows(Dy,Dx,Dt,win);
% optFlow = opticalFlowConstraint(data,V,win);
% 
% Nvx = size(V,3);
% Nvy = size(V,2);
% X = repmat(floor(win/2):win:Nvx*win,[Nvy,1]);
% Y = repmat((floor(win/2):win:Nvy*win)',[1,Nvx]);
% 
% figure;
% for t = 1:Nt
%     imagesc(squeeze(data(:,:,t)))
%     hold on
%     quiver(X,Y,squeeze(V(2,:,:,ceil(t/win))),squeeze(V(1,:,:,ceil(t/win))),'w')
%     hold off
%     pause()
% end
%% Test MATLAB optical flow
% clear opticFlow
% % opticFlow = opticalFlowLK('NoiseThreshold',0.0005);
% opticFlow = opticalFlowHS;
% h = figure;
% movegui(h);
% hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
% hPlot = axes(hViewPanel);
% 
% for t = 1:Nt
%     frameGray = squeeze(data(:,:,t));
%     flow = estimateFlow(opticFlow,frameGray);
%     imagesc(frameGray)
%     hold on
%     plot(flow,'DecimationFactor',[5 5],'ScaleFactor',1000,'Parent',hPlot);
%     q = findobj(gca,'type','Quiver');
%     q.Color = 'w';
%     hold off
%     pause()
% end

%% Compare my implementation of Horn-Schunk optical flow to MATLAB's
Fy = sobel(data,1);
Fx = sobel(data,2);
Ft = timeDiff(data);

u = zeros(size(data));
v = zeros(size(data));
avgKernel = [0 1 0; 1 0 1; 0 1 0];
smoothness = 1;

for i = 1:10
    % Update neighborhood averages
    uAvg = convn(u,avgKernel,'same');
    vAvg = convn(v,avgKernel,'same');

    u = uAvg - Fx.*(Fx.*uAvg + Fy.*vAvg + Ft)./(smoothness^2 + Fx.^2 + Fy.^2);
    v = vAvg - Fy.*(Fx.*uAvg + Fy.*vAvg + Ft)./(smoothness^2 + Fx.^2 + Fy.^2);
end

fig = figure;
fig.Position = [3.0326e+03 313.8000 1.6452e+03 554.8000];
gifDir = 'C:\Users\dpqb1\Documents\Outputs\opticalFlowFigs';
fName = 'opticalFlowExToy.gif';
opticFlow = opticalFlowHS;
for t = 1:80
    frameGray = squeeze(data(:,:,t));
    flow = estimateFlow(opticFlow,frameGray);
    subplot(2,1,1)
    imagesc(frameGray)
    hold on
    quiver(u(:,:,t),v(:,:,t))
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
    subplot(2,1,2)
    imagesc(frameGray)
    hold on
    quiver(flow.Vx,flow.Vy)
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off
    pause()
end
% sobelTrans(data,2)


