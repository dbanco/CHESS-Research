% Try different ways to compute the velocities in the coefficients
load('indepOutputs.mat')
X = gather(indepOutputs.X);
opt = indepOutputs.opt;
[N1,N2,KJ,T] = size(X);
K = 2;

plotOpticalFlow(X,K,opt)

Jterm = Ft == 0;
Jof = sum(vec(opticalFlowOp(Xtrue,u,v,K,0,Jterm)).^2);