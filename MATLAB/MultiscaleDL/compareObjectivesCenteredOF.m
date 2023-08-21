%% Load in different solutions and compute their objective values

% True solution
[y,y_true,K,J,N,M,T,Xtrue,Dtrue,scales] = gaus_linear_osc_signal_matched_small_zpad2_center(0.01);
center = round((M+1)/2);

[obj1,Jdf1,Jl11,Jof1,Jhs1] = computeObj(Dtrue,Xtrue,y,scales,center,N,M,K,5e-2,5e-1,1e-4);

% Indep solution
load(['C:\Users\dpqb1\Documents\Outputs\',...
     'toy3_center_exp_optFlow8_17_X0_D0_V0_zpad_HS0.002_sig_2\',...
     'output_j1_sig_1.00e-02_lam1_5.00e-02_lam2_0.00e+00.mat'])


[obj2,Jdf2,Jl12,Jof2,Jhs2] = computeObj(outputs.D,outputs.X,outputs.y,scales,center,N,M,K,5e-2,5e-1,1e-4);


% Coupled solution
load(['C:\Users\dpqb1\Documents\Outputs\',...
     'toy3_center_exp_optFlow8_17_X0_D0_V0_zpad_HS0.0001_sig_2\',...
     'output_j10_sig_1.00e-02_lam1_6.00e-02_lam2_5.00e-01.mat'])

[obj3,Jdf3,Jl13,Jof3,Jhs3] = computeObj(outputs.D,outputs.X,outputs.y,scales,center,N,M,K,5e-2,5e-1,1e-4);


%%
fprintf('Obj        Err       L1        OF        HS \n')
fprintf('%9.2e %9.2e %9.2e %9.2e %9.2e \n',obj1,Jdf1,Jl11,Jof1,Jhs1)
fprintf('%9.2e %9.2e %9.2e %9.2e %9.2e \n',obj2,Jdf2,Jl12,Jof2,Jhs2)
fprintf('%9.2e %9.2e %9.2e %9.2e %9.2e \n',obj3,Jdf3,Jl13,Jof3,Jhs3)