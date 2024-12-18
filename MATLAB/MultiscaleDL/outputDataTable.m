
% fNames = {'output_j1_1_sig_0.00e+00_lam1_1.00e-04_lam2_0.00e+00.mat',
%           'output_j5_1_sig_0.00e+00_lam1_1.00e-02_lam2_0.00e+00.mat',
%           'output_j10_1_sig_0.00e+00_lam1_5.33e-02_lam2_0.00e+00.mat',
%           'output_j100_1_sig_0.00e+00_lam1_7.53e-01_lam2_0.00e+00.mat'};

% Test 1.
% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_11_25_Dtrue1_Xtrue0\matched_log_results_sig_1\';

% Test 2.
% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_11_25_Dtrue1_Xzeros0\matched_log_results_sig_1\';

% Test 3.
% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_11_25_Dtrue0_Xzeros0\matched_log_results_sig_1\';

% Test 4.
topDir = 'C:\Users\dpqb1\Documents\Outputs2024_11_25_Dflat0_Xzeros0\matched_log_results_sig_1\';

% Get the list of all files in the directory
fileList = dir(topDir);

% Filter out directories and hidden files
isFile = ~[fileList.isdir]; % Logical array indicating files
fileList = fileList(isFile); % Keep only files

fileExtension = '.mat';
fNames = {fileList.name}'; % Get all filenames
fNames = fNames(endsWith(fNames, fileExtension));
numFiles = numel(fNames);

rowLabels = {'Lamba', 'Objective', 'Error', 'Penalty'};
numRows = numel(rowLabels);

% columnHeaders = {'Exp 1', 'Exp 2', 'Exp 3', 'Exp 4'};
table_vals = zeros(numRows,numFiles); % item, exp

for i = 1:numFiles
    load([topDir,fNames{i}])
    [Jfn,Jdf,Jl1,Jof,Jhs,lam_s] = computeObjMCDL(outputs);
    table_vals(1,i) = lam_s;
    table_vals(2,i) = Jfn;
    table_vals(3,i) = Jdf;
    table_vals(4,i) = Jl1;
end

% Sort table and to have increasing lambda values
[s_vals,s_inds] = sort(table_vals(1,:));
table_vals = table_vals(:,s_inds);

T = array2table(table_vals, 'RowNames', rowLabels);
disp(T);

writetable(T, 'tempTable.xlsx', 'WriteRowNames', true);
disp('Table exported to TableData.xlsx');