function sim_mcdl_reg_wrapper(lambdaVals,lambdaRegVals,j_s,j_reg,sigmas,i,opt,topDir,dataset)

figDir = [topDir,'_sig_',num2str(i)];
mkdir(figDir)

% Data  
if isfield(opt,'mcdl_file')
    [~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigmas(i),dataset);
    load(opt.mcdl_file)
    y = outputs.y;
else
    rng(1);
    % rng('shuffle');
    [y,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigmas(i),dataset);
    y = reshape(y,[1,N,T]);
end

% Model Setup
[K,J,scales] = datasetScales(dataset);
opt.DictFilterSizes = [1; M];
center = (M+1)/2;

lambda = lambdaVals(j_s);
lambda2 = lambdaRegVals(j_reg);

opt.lambda = lambda;
opt.lambda2 = lambda2;

% Initialize Coefficients and Dictionary
opt = initXD(opt,N,M,K,J,T,Xtrue,Dtrue);

if opt.mcdl_init > 0
    opt2 = opt;
    opt2.lambda2 = 0;
    opt2.MaxMainIter = opt.mcdl_init;
    [D,Y,~,~,~,~,~,~] = mcdl_regularized_gpu(opt2.G0,y,opt2,scales);
    opt.G0 = D;
    opt.Y0 = Y;
end

[D,Y,X,Dmin,Ymin,~,~,~] = mcdl_regularized_gpu(opt.G0,y,opt,scales);

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

AD = reSampleCustomArrayCenter3(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);

if opt.a_via_lam
    a = 0.95./outputs.lambda;
else
    a = opt.a;
end
metrics = compute_single_metrics(ADf,y,y_true,D,Dtrue,X,Xtrue,a,K,J,opt);
metrics.lambda = lambda;
metrics.lambda2 = lambda2;
outputs.metrics = metrics;

suffix = sprintf('_j%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e',...
                  j_s,j_reg,sigmas(i),outputs.lambda,outputs.lambda2);
save(fullfile(figDir,['output',suffix,'.mat']),'outputs');

% Generate figures
generateFiguresToy1zpad_center(figDir,outputs,suffix,[4,8],outputs.opt.useMin);

Bhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
plotDataRecon(y,Bhat,figDir,['y_recon',suffix,'.gif'])

close all

end