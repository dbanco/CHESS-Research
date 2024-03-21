function [Jfn,Jdf,Jl1,Jof,Jhs] = computeObj(outputs,scales,center,N,M,K,lam_s,lam_of,lam_hs)
%computeObj
D = outputs.D;
X = outputs.X;
S = outputs.y;
Uvel = outputs.Uvel;
Vvel = outputs.Vvel;

[AG,~] = reSampleCustomArrayCenter(N,D,scales,center);
AGpad = padarray(AG,[0 M-1 0 0],0,'post');
AGf = fft2(AGpad);
Xf = fft2(X);
% Jcn = norm(vec(Pcn(D) - D));
recon = unpad(ifft2(sum(bsxfun(@times,AGf,Xf),3),'symmetric'),M-1,'pre');
Jdf = sum(vec(abs(squeeze(recon)-squeeze(S)).^2))/2;
Jl1 = sum(abs(vec(X)));
% Jlg1 = rho*sum((X(:)-Y(:)+U(:)).^2);
% Jlg2 = sigma*sum((D(:)-G(:)+H(:)).^2);
if lam_of > 0
    if isempty(Uvel)
        [Uvel,Vvel,Fx,Fy,Ft] = computeHornSchunkDictPaperLS(X,K,[],[],lam_hs/lam_of,100);
    else
        [~,~,Fx,Fy,Ft] = computeHornSchunkDictPaperLS(X,K,[],[],lam_hs/lam_of,100);
    end
    [Jof, Jhs] = HSobjectivePaper(Fx,Fy,Ft,Uvel,Vvel,K,lam_hs);
else
    Jhs = 0;
    Jof = 0;
end
Jfn = Jdf + lam_s*Jl1 + lam_of*Jof + lam_hs*Jhs; % + Jlg1 + Jlg2;
Jfn = double(Jfn);


end