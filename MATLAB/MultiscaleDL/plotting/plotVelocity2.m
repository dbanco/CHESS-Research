K = 2;
[y,y_true,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched(0); %sigmas(1)

exps = {};

% alpha = 10, HSiters = 30, gamma = 0.1, no Jterm
exps{1} = ["C:\Users\dpqb1\Documents\Outputs\toy1_exp_optFlow2_noJterm_alpha10_sig_1\output_j5_sig_0.00e+00_lam1_2.00e-01_lam2_1.00e-01.mat"];
% alpha = 10, HSiters = 30, gamma = 1, no Jterm
exps{2} = ["C:\Users\dpqb1\Documents\Outputs\toy1_exp_optFlow2_noJterm_alpha10_sig_1\output_j6_sig_0.00e+00_lam1_2.00e-01_lam2_1.00e+00.mat"];
% alpha = 10, HSiters = 30, gamma = 5, no Jterm
exps{3} = ["C:\Users\dpqb1\Documents\Outputs\toy1_exp_optFlow2_noJterm_alpha10_sig_1\output_j1_sig_0.00e+00_lam1_2.00e-01_lam2_5.00e+00.mat"];
% alpha = 10, HSiters = 30, gamma = 10, no Jterm
exps{4} = ["C:\Users\dpqb1\Documents\Outputs\toy1_exp_optFlow2_noJterm_alpha10_sig_1\output_j2_sig_0.00e+00_lam1_2.00e-01_lam2_1.00e+01.mat"];
% alpha = 10, HSiters = 30, gamma = 20, no Jterm
exps{5} = ["C:\Users\dpqb1\Documents\Outputs\toy1_exp_optFlow2_noJterm_alpha10_sig_1\output_j3_sig_0.00e+00_lam1_2.00e-01_lam2_2.00e+01.mat"];
% alpha = 10, HSiters = 30, gamma = 100, no Jterm
exps{6} = ["C:\Users\dpqb1\Documents\Outputs\toy1_exp_optFlow2_noJterm_alpha10_sig_1\output_j4_sig_0.00e+00_lam1_2.00e-01_lam2_1.00e+02.mat"];

f = figure(3);
f.Position = [65 438 1.4136e+03 420];
cols = numel(exps) + 2;

load(exps{1})
opt = outputs.opt;
opt.Smoothness
Smoothness = opt.Smoothness;
load('indepOutputs.mat')
[uIndep,vIndep,~,~,~]  = computeHornSchunkDict(indepOutputs.X,K,Smoothness,opt.HSiters);

[u,v,Fy,Fx,Ft]  = computeHornSchunkDict(Xtrue,K,Smoothness,opt.HSiters);

% True solution
subplot(3,cols,1)
imagesc(squeeze(sum(v,3))')
title('True: v summed in time')
ylabel('scale')
xlabel('shift')

subplot(3,cols,2)
imagesc(squeeze(sum(vIndep,3))')
title('Indep:v summed in time')
ylabel('scale')
xlabel('shift')

subplot(3,cols,cols+1)
imagesc(squeeze(sum(u,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')

% Indep Solution
subplot(3,cols,cols+2)
imagesc(squeeze(sum(uIndep,3))')
title('u summed in time')
ylabel('scale')
xlabel('shift')

subplot(3,cols,2*cols+1)
imagesc(squeeze(sum(Xtrue,2)))
title('Atom usage in time')
ylabel('scale')
xlabel('time')  

subplot(3,cols,2*cols+2)
imagesc(squeeze(sum(indepOutputs.X,2)))
title('Atom usage in time')
ylabel('scale')
xlabel('time')

% Coupled Solutions
for i = 1:numel(exps)
    load(exps{i})
    if ~isfield(outputs,'Uvel') 
        [Uvel,Vvel,~,~,~]  = computeHornSchunkDict(outputs.X,K,opt.Smoothness,opt.HSiters);
        outputs.Uvel = Uvel;
        outputs.Vvel = Vvel;
    end
    
    subplot(3,cols,2+i)
    imagesc(squeeze(sum(outputs.Vvel,3))')
    title('Coupled: v summed in time')
    ylabel('scale')
    xlabel('shift')
    
    subplot(3,cols,cols+2+i)
    imagesc(squeeze(sum(outputs.Uvel,3))')
    title('u summed in time')
    ylabel('scale')
    xlabel('shift')
      
    subplot(3,cols,2*cols+2+i)
    imagesc(squeeze(sum(outputs.X,2)))
    title('Atom usage in time')
    ylabel('scale')
    xlabel('time')
end




