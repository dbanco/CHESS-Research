baseFileName = 'coupled_fit_%i_%i.mat';
dir1 = 'E:\MMPAD_data\ring1_zero_coupled_CG_TVphi4';
dir2 = 'E:\MMPAD_data\ring1_zero_coupled_CG_TVphi4b';
outdir = 'E:\MMPAD_data\ring1_zero_coupled_CG_TVphi4c';
mkdir(outdir)

% Define new set of lambda_values
load(fullfile(dir1,sprintf(baseFileName,1,1)))
lambda_values = P.lambda_values;
num_ims = P.num_ims;
M1 = numel(P.lambda_values);

load(fullfile(dir2,sprintf(baseFileName,1,1)))
lambda_values = [lambda_values P.lambda_values];
M2 = numel(P.lambda_values);

% Move all files
for k = 1:num_ims
    j = 1;
    % Transfer dir1
    for i = 1:M1
        f1 = fullfile(dir1,sprintf(baseFileName,i,k));
        f2 = fullfile(outdir,sprintf(baseFileName,j,k));
        movefile(f1,f2)
        j = j + 1;
    end
    % Transfer dir2
    for i = 1:M2
        f1 = fullfile(dir2,sprintf(baseFileName,i,k));
        f2 = fullfile(outdir,sprintf(baseFileName,j,k));
        movefile(f1,f2)
        j = j + 1;
    end
end

%% Correct lambda values
lambda_values = [logspace(-6,-4.1,15) logspace(-4,1,30)];
P.lambda_values = lambda_values;
for k = 1:num_ims
    for i = 1:M1
        load(fullfile(outdir,sprintf(baseFileName,i,k)))
        P.lambda_values = lambda_values;
        save(fullfile(outdir,sprintf(baseFileName,i,k)),'err','l1_norm','obj','P','x_hat')
    end
end