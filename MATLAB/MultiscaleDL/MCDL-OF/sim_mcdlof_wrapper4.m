function sim_mcdlof_wrapper4(lambdaVals,lambdaOFVals,lambdaHSVals,j_s,j_of,j_hs,sigmas,i,opt,topDir,dataset,K,scales)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Optical flow coupled solution

figDir = [topDir,'_sig_',num2str(i)];
mkdir(figDir)

% Data  
if isfield(opt,'mcdl_file')
    [~,~,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigmas(i),dataset);
    load(opt.mcdl_file)
    y = outputs.y;
else
    rng('shuffle');
    [y,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigmas(i),dataset);
    y = reshape(y,[1,N,T]);
end

% Model Setup
J = size(scales{1},2);
opt.DictFilterSizes = [1; M];
center = (M+1)/2;

lambda = lambdaVals(j_s);
lambda2 = lambdaOFVals(j_of);
opt.Smoothness = lambdaHSVals(j_hs);

% Initialize Coefficients and Dictionary
opt = initXD(opt,N,M,K,J,T,Xtrue,Dtrue);

if opt.mcdl_init
    opt2 = opt;
    [Uvel,Vvel,~,~,~] = computeHornSchunkDictPaperLS(opt2.Y0,K,[],[],opt2.Smoothness/lambda2,opt2.HSiters);
    opt2.UpdateVelocity = 1;
    opt2.Smoothness = 0;
    opt2.MaxMainIter = 500;
    [D,Y,X,Dmin,Ymin,~,~,~,~,~] = cbpdndl_cg_OF_multiScales_gpu_zpad_center3(opt2.G0, y, lambda,lambda2, opt2, scales,Uvel,Vvel);
    opt.G0 = Dmin;
end

[Uvel,Vvel,~,~,~] = computeHornSchunkDictPaperLS(opt.Y0,K,[],[],opt.Smoothness/lambda2,opt.HSiters);
opt.UpdateVelocity = 1;
[D,Y,X,Dmin,Ymin,Uvel,Vvel,~,~,~] = cbpdndl_cg_OF_multiScales_gpu_zpad_center3(opt.G0, y, lambda,lambda2, opt, scales,Uvel,Vvel);

% Save outputs
outputs = struct();
outputs.y = y;
outputs.D = D;
outputs.X = X;
outputs.Y = Y;
outputs.Dmin = Dmin;
outputs.Ymin = Ymin;
outputs.scales = scales;
outputs.N = N;
outputs.M = M;
outputs.T = T;
outputs.K = K;
outputs.J = J;
outputs.opt = opt;
outputs.lambda = lambda;
outputs.lambda2 = lambda2;
outputs.lambda3 = opt.Smoothness;
outputs.Uvel = Uvel;
outputs.Vvel = Vvel;
suffix = sprintf('_j%i_%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e_lam3_%0.2e',...
                  j_s,j_of,j_hs,sigmas(i),outputs.lambda,outputs.lambda2,outputs.lambda3);
save(fullfile(figDir,['output',suffix,'.mat']),'outputs');

% Generate figures
generateFiguresToy1zpad_center(figDir,outputs,suffix,[4,8],outputs.opt.useMin);

AD = reSampleCustomArrayCenter3(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Bhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
plotDataRecon(y,Bhat,figDir,['y_recon',suffix,'.gif'])

close all

end