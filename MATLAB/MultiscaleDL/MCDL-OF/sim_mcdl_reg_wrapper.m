function sim_mcdl_reg_wrapper(lambdaVals,lambdaOFVals,j_s,j_of,sigmas,i,opt,topDir,dataset,K,scales)
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

% Initialize Coefficients and Dictionary
opt = initXD(opt,N,M,K,J,T,Xtrue,Dtrue);

if opt.mcdl_init
    opt2 = opt;
    opt2.MaxMainIter = 500;
    [D,Y,X,Dmin,Ymin,~,~,~] = mcdl_regularized_gpu(opt2.G0,y,lambda,0,opt2,scales);
    opt.G0 = Dmin;
end

[D,Y,X,Dmin,Ymin,~,~,~] = mcdl_regularized_gpu(opt.G0,y,lambda,lambda2,opt,scales);

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
suffix = sprintf('_j%i_%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e',...
                  j_s,j_of,j_hs,sigmas(i),outputs.lambda,outputs.lambda2);
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