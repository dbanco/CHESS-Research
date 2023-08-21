% Load true solution and indep solution
[y,y_true,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched(0);
load('indepOutputs.mat')
Xindep = gather(indepOutputs.X);
K = 2;

% Compute current optical flow objective with different parameters
% tic
% [u1,v1,~,~,~] = computeHornSchunkDict(Xtrue,K,10,30);
% [u2,v2,~,~,~] = computeHornSchunkDict(Xindep,K,10,30);
% toc
% 
% figure(11)
% subplot(2,1,1)
% imagesc(squeeze(sum(v2,3))')
% title('v summed in time')
% ylabel('scale')
% xlabel('shift')
% 
% subplot(2,1,2)
% imagesc(squeeze(sum(u2,3))')
% title('u summed in time')
% ylabel('scale')
% xlabel('shift')

% Try other options for optical flow computation and parameters
% OFtrue1 = norm(opticalFlowOp(Xtrue,u1,v1,K),'fro');
% OFindep1 = norm(opticalFlowOp(Xindep,u2,v2,K),'fro');

%% Velocities frame-by-frame
opt.Smoothness = 1e-4;
smth = opt.Smoothness;
opt.HSiters = 300;
maxIt = opt.HSiters;
% plotOpticalFlow3(Xtrue,K,opt)

%% Different optical flows
% [u1,v1,Fx1,Fy1,Ft1] = computeHornSchunkDict(Xtrue,K,smth,maxIt);
% [u2,v2,Fx2,Fy2,Ft2] = computeHornSchunkDictLS(Xtrue,K,[],[],smth,maxIt);
[u3,v3,Fx3,Fy3,Ft3] = computeHornSchunkDictPaper(Xtrue,K,smth,maxIt);
[u4,v4,Fx4,Fy4,Ft4] = computeHornSchunkDictPaperLS(Xtrue,K,[],[],smth,maxIt);

% OFtrue1 =     norm(opticalFlowOp(Xtrue,u1,v1,K),'fro');
% OFtrueLS1 =   norm(opticalFlowOp(Xtrue,u2,v2,K),'fro');
% OftruePaper = norm(opticalFlowOpPaper(Xtrue,u3,v3,K),'fro');
OftruePaperLS = norm(opticalFlowOpPaper(Xtrue,u4,v4,K),'fro');

% [u1i,v1i,~,~,~] = computeHornSchunkDict(Xindep,K,smth,maxIt);
% [u2i,v2i,~,~,~] = computeHornSchunkDictLS(Xindep,K,[],[],smth,maxIt);
% [u3i,v3i,Fx3i,Fy3i,Ft3i] = computeHornSchunkDictPaper(Xindep,K,smth,maxIt);
[u4i,v4i,Fx4i,Fy4i,Ft4i] = computeHornSchunkDictPaperLS(Xindep,K,[],[],smth,maxIt);

% OFindep1 =     norm(opticalFlowOp(Xindep,u1i,v1i,K),'fro');
% OFindepLS1 =   norm(opticalFlowOp(Xindep,u2i,v2i,K),'fro');
% OfindepPaper = norm(opticalFlowOpPaper(Xindep,u3i,v3i,K),'fro');
OfindepPaperLS = norm(opticalFlowOpPaper(Xindep,u4i,v4i,K),'fro');

% [objOF3, objHS3, sys_3] = HSobjectivePaper(Fx3,Fy3,Ft3,u3,v3,K,smth);
[objOF4, objHS4, sys_4] = HSobjectivePaper(Fx4,Fy4,Ft4,u4,v4,K,smth);
[objOF4i, objHS4i, sys_4i] = HSobjectivePaper(Fx4i,Fy4i,Ft4i,u4i,v4i,K,smth);

ofTrue = objOF4 + objHS4*smth;
ofIndep = objOF4i + objHS4i*smth;

% sys_3norm = norm(sys_3(:));
sys_4norm = norm(sys_4(:));
sys_4inorm = norm(sys_4i(:));

figure(4)
subplot(2,1,1)
imagesc(squeeze(sum(v4,3))')
title('v summed in time')
ylabel('scale')
xlabel('shift')

subplot(2,1,2)
imagesc(squeeze(sum(u4,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')

figure(41)
subplot(2,1,1)
imagesc(squeeze(sum(v4i,3))')
title('v summed in time')
ylabel('scale')
xlabel('shift')

subplot(2,1,2)
imagesc(squeeze(sum(u4i,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')


%% Inspect constraint equation terms
t = 7;
window1 = 140:170;
window2 = 1:20;

fig1 = inspectOFconstraint(u1,v1,Fx1,Fy1,Ft1,window1,window2,t);
fig2 = inspectOFconstraint(u2,v2,Fx2,Fy2,Ft2,window1,window2,t);
fig3 = inspectOFconstraint(u3,v3,Fx3,Fy3,Ft3,window1,window2,t);
fig4 = inspectOFconstraint(u4,v4,Fx4,Fy4,Ft4,window1,window2,t);
fig4i = inspectOFconstraint(u4i,v4i,Fx4i,Fy4i,Ft4i,window1,window2,t);



%% Plot summed velocities
figure(1)
subplot(2,1,1)
imagesc(squeeze(sum(v1,3))')
title('v summed in time')
ylabel('scale')
xlabel('shift')

subplot(2,1,2)
imagesc(squeeze(sum(u1,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')

figure(2)
subplot(2,1,1)
imagesc(squeeze(sum(v2,3))')
title('v summed in time')
ylabel('scale')
xlabel('shift')

subplot(2,1,2)
imagesc(squeeze(sum(u2,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')

figure(3)
subplot(2,1,1)
imagesc(squeeze(sum(v3,3))')
title('v summed in time')
ylabel('scale')
xlabel('shift')

subplot(2,1,2)
imagesc(squeeze(sum(u3,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')

figure(4)
subplot(2,1,1)
imagesc(squeeze(sum(v4,3))')
title('v summed in time')
ylabel('scale')
xlabel('shift')

subplot(2,1,2)
imagesc(squeeze(sum(u4,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')

figure(41)
subplot(2,1,1)
imagesc(squeeze(sum(v4i,3))')
title('v summed in time')
ylabel('scale')
xlabel('shift')

subplot(2,1,2)
imagesc(squeeze(sum(u4i,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')


% Compare sobel and sobel transposed
% d1 = sobel(squeeze(Xindep),1);
% d2 = sobel(squeeze(Xindep),1,1);
% norm(d1(:)+d2(:))
% figure
% subplot(2,1,1)
% imagesc(d1(:,:,2))
% subplot(2,1,2)
% imagesc(d2(:,:,2))


