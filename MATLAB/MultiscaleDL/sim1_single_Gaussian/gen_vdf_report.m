% Ensure the Report Generator Toolbox is installed.

% Create a Word document
import mlreportgen.dom.*;
doc = Document('Vdf_report', 'docx');
open(doc);

dataset = 'steps_matched';

testInds = [1,3,5,10,50,90];
lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];


dfixes = {1,0};
xinits = {'true','zeros'};
dinits = {'true','flat'};
penalties = {'l1-norm','log'};

for sig = [1,3]
    for s3 = 1:2
        for s4 = 1:2
            for s2 = 1:2
                for s1 = 1:2
                    
if (dfixes{s3} == 1) && strcmp(dinits{s2}, 'flat')
    continue
end
% Determine directory
topDir = ['E:\Outputs2024_12_3',...
                    '_D',dinits{s2},num2str(dfixes{s3}),...
                    '_X',xinits{s4},num2str(0),'\',...
                    dataset,'_',penalties{s1},'_results_sig_',num2str(sig)];

% Add a title
if dfixes{s3}
    fixStr = ' (fixed)';
else
    fixStr = '';
end
paraStr = [dataset,', ',penalties{s1},...
            ', initialized D as ',dinits{s2},fixStr,...
            ', initialized X as ',xinits{s4}];
title = Paragraph(paraStr);
title.Style = {Bold(true), FontSize('14pt')};
append(doc, title);

% Get each image
% Create a table for 2x3 layout
imageTable = Table();
imageTable.Width = '100%'; % Adjust width as needed

% Style for the table
rowStyle = {RowHeight('1.5in')}; % Adjust the height of each row
cellStyle = {HAlign('center'), VAlign('middle')};

% Setup for data table
rowLabels = {'Lamba', 'Objective', 'Error', 'Penalty'};
numRows = numel(rowLabels);
table_vals = zeros(numRows,length(testInds));

% Iterate over testInds and add images
row = TableRow(); % Start a new row
for idx = 1:length(testInds)
    i = testInds(idx);
    try
        imageFile = dir(fullfile(topDir, ['vdf_j', num2str(i), '_*.png'])); % Adjust extension as needed
        imPath = fullfile(imageFile.folder, imageFile.name);
        imInfo = imfinfo(imPath);
        
        % Create and resize image
        originalWidth = imInfo.Width;
        originalHeight = imInfo.Height;
        aspectRatio = originalWidth / originalHeight;
        desiredHeightInches = 1; % 1 inch
        desiredWidthInches = desiredHeightInches * aspectRatio;

        image = Image(imPath);
        image.Height = sprintf('%fin', desiredHeightInches);
        image.Width = sprintf('%fin', desiredWidthInches);
        
        % Create a table cell and add the image
        cell = TableEntry();
        cell.Style = cellStyle; % Center alignment
        append(cell, image);

        % Add caption
        lambdaSymbol = char(955); % Unicode for Greek lowercase lambda (λ)
        caption = Paragraph(sprintf('%s = %0.2e',lambdaSymbol, lambdaVals(i)));
        caption.Style = {Italic(true), FontSize('10pt')};
        append(cell, caption);

        append(row, cell); % Add cell to the row

        % Check if the row is complete (3 cells)
        if mod(idx, 3) == 0
            append(imageTable, row); % Add the completed row to the table
            row = TableRow(); % Start a new row
        end

        % Get table data
        dataFile = dir(fullfile(topDir, ['output_j', num2str(i), '_*.mat']));
        dataPath = fullfile(dataFile.folder, dataFile.name);
        load(dataPath)

        [Jfn,Jdf,Jl1,Jof,Jhs,lam_s] = computeObjMCDL(outputs);
        table_vals(1,idx) = lam_s;
        table_vals(2,idx) = Jfn;
        table_vals(3,idx) = Jdf;
        table_vals(4,idx) = Jl1;

    catch
        fprintf('File not found for index %d\n', i);
    end
end

% If there's an incomplete row, add it to the table
if ~isempty(row.Children)
    append(imageTable, row);
end

% Append the table to the document
append(doc, imageTable);

% Add a line break
append(doc, Paragraph(''));

 % Sort table and to have increasing lambda values
[s_vals,s_inds] = sort(table_vals(1,:));
table_vals = table_vals(:,s_inds);

T = array2table(table_vals, 'RowNames', rowLabels);% Create a DOM table
domTable = Table();

% Add rows with row labels
for rowIdx = 1:height(T)
    row = TableRow();
    % Add row label
    append(row, TableEntry(Paragraph(T.Properties.RowNames{rowIdx})));
    % Add table data
    for colIdx = 1:width(T)
        value = T{rowIdx, colIdx};
        entry = TableEntry(Paragraph(num2str(value, '%.3f'))); % Format values to 3 decimal places
        append(row, entry);
    end
    append(domTable, row);
end

% Style the table
domTable.Border = 'solid';
domTable.BorderWidth = '1pt';
domTable.TableEntriesHAlign = 'left';

% Add table to document
append(doc, domTable);

                end
            end
        end
    end
end

% Close the document
close(doc);

