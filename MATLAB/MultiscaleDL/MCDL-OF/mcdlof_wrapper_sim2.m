function mcdlof_wrapper_sim2(lambdaVals,lambdaOFVals,lambdaHSVals,j_s,j_of,j_hs,sigmas,i,opt,K,scales,topDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Optical flow coupled solution

figDir = [topDir,'_sig_',num2str(i)];
mkdir(figDir)

% Data  
[y,~,N,M,T] = generate_signal_pair(sigmas(i));
y = reshape(y,[1,N,T]);
center = (M+1)/2;

lambda = lambdaVals(j_s);
lambda2 = lambdaOFVals(j_of);
opt.Smoothness = lambdaHSVals(j_hs);

[Uvel,Vvel,~,~,~] = computeHornSchunkDictPaperLS(opt.Y0,K,[],[],opt.Smoothness/lambda2,opt.HSiters);
opt.UpdateVelocity = 1;
[D,X,Dmin,Xmin,Uvel,Vvel,~,~,~] = cbpdndl_cg_OF_multiScales_gpu_zpad_center(opt.G0, y, lambda,lambda2, opt, scales,Uvel,Vvel);

% Save outputs
outputs = struct();
outputs.y = y;
outputs.D = D;
outputs.X = X;
outputs.Dmin = Dmin;
outputs.Xmin = Xmin;
outputs.scales = scales;
outputs.N = N;
outputs.M = M;
outputs.T = T;
outputs.K = K;
outputs.opt = opt;
outputs.lambda = lambda;
outputs.lambda2 = lambda2;
outputs.Uvel = Uvel;
outputs.Vvel = Vvel;
suffix = sprintf('_j%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e',...
                  j_of,sigmas(i),outputs.lambda,outputs.lambda2);
save(fullfile(figDir,['output',suffix,'.mat']),'outputs');

% Generate figures
generateFiguresToy1zpad_center(figDir,outputs,suffix,[4,8]);
%         generateFiguresToy1min([figDir,'min'],outputs,suffix)
%         generateFiguresToy1([figDir,'indep'],inde,suffix)

AD = reSampleCustomArrayCenter(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
plotDataRecon(y,Yhat,figDir,['y_recon',suffix,'.gif'])
close all

end