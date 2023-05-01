% Load outputs
outputDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_200exact_sig_1';
K = 2;
smoothness = 1e-8;
maxIt = 1;
close all

% True solution
load('C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_200exact_sig_1\output_j1_sig_0.00e+00_lam1_0.00e+00_lam2_0.00e+00.mat');
data = squeeze(outputs.X);
[N, KJ, T] = size(data);
[u,v,Fy,Fx,Ft]  = computeHornSchunkDict(data,K,smoothness,maxIt);
term = data;
term(Ft==0) = 0;
of = abs(Fx.*u + Fy.*v + Ft);
sum( abs(Fx.*u + Fy.*v + Ft),'all')
sum( abs(Fx.*u + Fy.*v + Ft + term),'all')
sum( abs(Fx.*u + Fy.*v + Ft+ term)>1e-8,'all')


% figure
% plot(squeeze(sum(sum(of,1),2)))
% velCon(data,K,u,v)

% L1-norm solution
load('C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_initD0_200_sig_1\output_j1_sig_0.00e+00_lam1_1.50e-01_lam2_0.00e+00.mat');
data = squeeze(outputs.X);
sum( abs(opticalFlowOp(data,u,v,K)),'all')
sum( abs(opticalFlowOp(data,u,v,K,0,1)),'all')
sum( abs(opticalFlowOp(data,u,v,K))>0,'all')


[u,v,Fy,Fx,Ft]  = computeHornSchunkDict(data,K,smoothness,maxIt);
sum( abs(opticalFlowOp(data,u,v,K)),'all')
sum( abs(opticalFlowOp(data,u,v,K,0,1)),'all')
sum( abs(opticalFlowOp(data,u,v,K))>0,'all')


%% Compute opticalflow
[u,v,Fy,Fx,Ft]  = computeHornSchunkDict(data,K,smoothness,maxIt);
OF = sum( abs(opticalFlowOp(data,u,v,K)),'all');
% u = ones(size(u));
% v = ones(size(v));

% View velocities in scale/shift
figure
subplot(1,2,1)
imagesc(squeeze(sum(u,3)))
title('u (scale) summed in time')
ylabel('shift')
xlabel('scale')

subplot(1,2,2)
imagesc(squeeze(sum(v,3)))
title('v (shift) summed in time')
ylabel('shift')
xlabel('scale')

figure
imagesc(squeeze(sum(data,3)))
title('coefficients summed in time')
ylabel('shift')
xlabel('scale')

%% 
winY = 155:160;
winX = 12:18;
figure
for t = 9
    subplot(5,1,1)
    imagesc(Fx(winY,winX,t)')
    title('F_x')
    subplot(5,1,2)
    imagesc(u(winY,winX,t)')
    title('u')
    subplot(5,1,3)
    imagesc(Fy(winY,winX,t)')
    title('F_y')
    subplot(5,1,4)
    imagesc(v(winY,winX,t)')
    title('v')
    subplot(5,1,5)
    imagesc(Ft(winY,winX,t)')
    title('F_t')
end

figure
for t = 9
    subplot(5,1,1)
    imagesc(Fx(winY,winX,t)'.*u(winY,winX,t)')
    title('u*F_x')
    subplot(5,1,2)
    imagesc(Fy(winY,winX,t)'.*v(winY,winX,t)')
    title('v*F_y')
    subplot(5,1,3)
    imagesc(Ft(winY,winX,t)'+Fy(winY,winX,t)'.*v(winY,winX,t)'+Fx(winY,winX,t)'.*u(winY,winX,t)')
    title('OF')
    subplot(5,1,4)
    imagesc(data(winY,winX,t)')
    title('x_{t}')
    subplot(5,1,5)
    imagesc(data(winY,winX,t-1)')
    title('x_{t+1}')
end

%%
figure
imagesc(squeeze(sum(u.*Fx + v.*Fy,3)))

figure
subplot(1,2,1)
imagesc(squeeze(sum(Fy.*v,3)))
title('v*F_y summed in time')
ylabel('shift')
xlabel('scale')

subplot(1,2,2)
imagesc(squeeze(sum(Fx.*u,3)))
title('u*F_x summed in time')
ylabel('shift')
xlabel('scale')



figure
subplot(1,2,1)
imagesc(squeeze(sum(Fy,3)))
title('F_y summed in time')
ylabel('shift')
xlabel('scale')

subplot(1,2,2)
imagesc(squeeze(sum(Fx,3)))
title('F_x summed in time')
ylabel('shift')
xlabel('scale')

figure
subplot(1,2,1)
imagesc(squeeze(sum(Ft,1)))
title('F_t summed in space')
ylabel('scale')
xlabel('time')

subplot(1,2,2)
imagesc(squeeze(sum(Ft,2)))
title('F_t summed in scale')
ylabel('shift')
xlabel('time')

%% View optical flow
fig = figure;
fig2 = figure;
fig.Position = [2716.2 313.8 1.6452e+03 554.8000];
% outputDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_200exact_sig_1';
outputDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_initD0_200_sig_1'
gifDir = outputDir;
fName = 'opticalFlow_1min.gif';
opticFlow = opticalFlowHS;
for t = 1:T
    frameGray = squeeze(data(:,:,t));
%     flow = estimateFlow(opticFlow,frameGray');
    figure(fig)
    imagesc(frameGray')
    hold on
%     plot(flow,'DecimationFactor',[5 5],'ScaleFactor',60,'Parent',hPlot);
%     quiver(flow.Vx,flow.Vy)
    quiver(v(:,:,t)',u(:,:,t)')
    q = findobj(gca,'type','Quiver');
    q.Color = 'w';
    hold off
    title(sprintf('t=%i',t))
    
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    figure(fig2)
    imagesc(of(:,:,t)')

%     % Save the image to a GIF file
%     if t == 1
%         outFile = fullfile(gifDir,fName);
%         imwrite(imind,cm,outFile,'gif','Loopcount',inf,'DelayTime',0.2);
%     else
%         imwrite(imind,cm,outFile,'gif','WriteMode','append','DelayTime',0.2);
%     end

    pause(0.0001)
end