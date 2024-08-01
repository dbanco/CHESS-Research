function outputs = saveOutputsMS(y,D,X,Dmin,Xmin,Uvel,Vvel,lambda,lambda2,lambda3,N,M,K,T,scales,center,opt,j_s,j_of,j_hs,outDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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
            outputs.lambda3 = lambda3;
            outputs.Uvel = Uvel;
            outputs.Vvel = Vvel;
            suffix = sprintf('_j_%i_%i_%i_lam1_%0.2e_lam2_%0.2e_lam3_%0.2e',...
                              j_s,j_of,j_hs,lambda,lambda2,lambda3);
            save(fullfile(outDir,['output',suffix,'.mat']),'outputs');

            % Generate figures
            generateFiguresToy1zpad_center(outDir,outputs,suffix,[12,4])

            AD = reSampleCustomArrayCenter(N,D,scales,center);
            AD = padarray(AD,[0 M-1 0],0,'post');
            ADf = fft2(AD);
            Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
            plotDataRecon(y,Yhat,outDir,['y_recon',suffix,'.gif'])
            close all
end