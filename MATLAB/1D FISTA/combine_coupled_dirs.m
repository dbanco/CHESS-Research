baseFileName = 'coupled_fit_%i.mat';
dir1 = 'E:\MMPAD_data\ring1_zero_coupled_CG_TVphi4b';
dir2 = 'E:\MMPAD_data\ring1_zero_coupled_CG_TVphi4';
outdir = 'E:\MMPAD_data\ring1_zero_coupled_CG_TVphi4c';
mkdir(outdir)

% Define new set of lambda_values
load(fullfile(dir1,sprintf(baseFileName,11)))
lambda_values = P.lambda_values;
num_ims = P.num_ims;
M1 = numel(P.lambda_values);
M1 = 15
load(fullfile(dir2,sprintf(baseFileName,11)))
lambda_values = [lambda_values P.lambda_values];
M2 = numel(P.lambda_values);
M2 = 30
j = 1;
% Transfer dir1
for i = 1:M1
    f1 = fullfile(dir1,sprintf(baseFileName,i));
    f2 = fullfile(outdir,sprintf(baseFileName,j));
    movefile(f1,f2)
    j = j + 1;
end
% Transfer dir2
for i = 1:M2
    f1 = fullfile(dir2,sprintf(baseFileName,i));
    f2 = fullfile(outdir,sprintf(baseFileName,j));
    movefile(f1,f2)
    j = j + 1;
end

%% Correct lambda values
lambda_values = [logspace(-6,-4.1,15) logspace(-4,1,30)];
P.lambda_values = lambda_values;
for i = 1:M1
    load(fullfile(outdir,sprintf(baseFileName,i)))
    P.lambda_values = lambda_values;
    save(fullfile(outdir,sprintf(baseFileName,i)),'err','l1_norm','obj','P','x_hat')
end