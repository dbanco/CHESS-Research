function hexd_realTime_gui(filePath,outPath,r,r1,r2,center,omInd,omInd2,startFnum)
    % Create the GUI figure
    fig = uifigure('Name', 'hexd_realTime_gui', 'Position', [10 40 1500 520]);
    
    % Create UI components
    startButton = uibutton(fig, 'Text', 'Start/Resume', 'Position', [20 20 120 40]);
    
    pauseButton = uibutton(fig, 'Text', 'Pause', 'Position', [160 20 120 40]);
    continueBox = uicheckbox(fig, 'Text', 'Continue','Value',1,'Position',[285 25 20 20]);
    resetButton = uibutton(fig, 'Text', 'Reset', 'Position', [300 20 120 40]);
    processH5Button = uibutton(fig, 'Text', 'Read H5', 'Position', [500 20 120 40]);
    nextImageButton = uibutton(fig, 'Text', 'Next Image', 'Position', [640 20 120 40]);
    loadGE2Button = uibutton(fig, 'Text', 'Read Ge2', 'Position', [780 20 120 40]);
        
    fileLabel = uilabel(fig, 'Text', 'No file being processed', 'Position', [20 495 150 30]);
    

    lambdaField = uieditfield(fig,'numeric','Position',[190 495 60 30],'Limits',[0 100]);
    rad1Field = uieditfield(fig,'numeric','Position',[295 495 40 30],'Limits',[0 10000]);
    rad2Field = uieditfield(fig,'numeric','Position',[345 495 40 30],'Limits',[0 10000]);
    center1Field = uieditfield(fig,'numeric','Position',[475 495 40 30],'Limits',[0 10000]);
    center2Field = uieditfield(fig,'numeric','Position',[525 495 40 30],'Limits',[0 10000]);
    updateButton = uibutton(fig,'Text','Update','Position',[600 495 60 30]);

    uilabel(fig, 'Text', 'lambda:', 'Position', [170 495 40 30]);
    uilabel(fig, 'Text', 'radii:','Position', [260 495 50 30]);
    uilabel(fig, 'Text', 'center:', 'Position', [400 495 70 30]);
    
    % Set up axes
    dictionaryAxes = uiaxes(fig, 'Position', [20 300 600 200]);
    vdfAxes = uiaxes(fig, 'Position', [620 300 400 200]);
    awmvPlot = uiaxes(fig, 'Position', [1020 300 400 200]);
    bImageAxes = uiaxes(fig, 'Position', [20 180 1400 120]);
    bhatImageAxes = uiaxes(fig, 'Position', [20 60 1400 120]);

    axis(dictionaryAxes, 'tight');
    title(dictionaryAxes,'Dictionary')
    xlabel(dictionaryAxes, '\eta');
    ylabel(dictionaryAxes, '2\theta');
    set(dictionaryAxes,'xtick',[],'ytick',[])

    axis(vdfAxes, 'tight');
    title(vdfAxes,'VDF')
    xlabel(vdfAxes, '\eta');
    ylabel(vdfAxes, '2\theta');
    set(vdfAxes,'xtick',[],'ytick',[])

    title(awmvPlot,'\eta and 2\theta AWMV(t)')
    xlabel(awmvPlot, 'Image number');
    ylabel(awmvPlot, 'AWMV');

    axis(bImageAxes, 'tight');
    title(bImageAxes,'Data')
    xlabel(bImageAxes, '\eta');
    ylabel(bImageAxes, '2\theta');
    set(bImageAxes,'xtick',[],'ytick',[])

    axis(bhatImageAxes, 'tight');
    title(bhatImageAxes,'Reconstruction')
    xlabel(bhatImageAxes, '\eta');
    ylabel(bhatImageAxes, '2\theta');
    set(bhatImageAxes,'xtick',[],'ytick',[])
        
    fileName = [];
    processedFiles = {};
    printDot = 0;

        % GE2 detector/ring position
%     center = [1025,1020];
%     r1 = 430;
%     r2 = 445;

    % Test Dexela detector/ring position
%     filePath = '/nfs/chess/raw/2022-3/id1a3/miller-3528-a/ff-c103-90-s2-1/2/ff';
%     r1 = 1195;
%     r2 = 1225;
%     center = [1944, 3560];
    
    % Dexela detector/ring position

    filePattern = fullfile(filePath, '/*/ff/ff1*.h5'); 
    filePattern2 = fullfile(filePath, '/*/ff/ff2*.h5'); 
    mkdir(outPath)
    
    % Parameters
    % Basis function variance parameters
    P.K1 = 10;
    P.K2 = 15;
    P.K = P.K1*P.K2;
    P.sigma1 = linspace(0.5,  4,    P.K1);
    P.sigma2 = linspace(0.5,  4,   P.K2); % was 12
    P.basis = 'norm2';
 
    
    % fista params
    params.lambda = 3e-3; % sparsity penalty
    lambdaField.Value = params.lambda;
    params.L = 200;  %
    params.beta = 2; %
    params.noBacktrack = 0;
    
    params.isNonnegative = 1; % flag to enforce nonnegativity
    
    params.stoppingCriterion = 'OBJECTIVE_VALUE';
    params.maxIter = 1000;
    params.tolerance = 1e-2;
    
    params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
    params.verbose = 1;      % flag to print objective values at each iteration 
    P.params = params;

    rad1Field.Value = r1;
    rad2Field.Value = r2;
    center1Field.Value = center(1);
    center2Field.Value = center(2);

    x_init = [];
    D = []; 
    Df = [];
    b = [];
    normb = 1;
    vdf = [];
    awmv_eta = [];
    awmv_2th = [];
    t = 0;  % Number of files processed
    rel_error = 1;
    oldbsize = [];
    
    
    % Set up timer
    timerObj = timer('ExecutionMode', 'fixedRate', 'Period', 1, 'TimerFcn', @timerCallback);
    
    % Button callbacks
    startButton.ButtonPushedFcn = @(~,~) startTimer();
    pauseButton.ButtonPushedFcn = @(~,~) pauseTimer();
    resetButton.ButtonPushedFcn = @(~,~) resetData();
    processH5Button.ButtonPushedFcn = @(~,~) processH5();
    nextImageButton.ButtonPushedFcn = @(~,~) nextImage();
    loadGE2Button.ButtonPushedFcn = @(~,~) processGE2();
    lambdaField.ValueChangedFcn = @(~,~) updateFields();
    rad1Field.ValueChangedFcn = @(~,~) updateFields();
    rad2Field.ValueChangedFcn = @(~,~) updateFields();
    center1Field.ValueChangedFcn = @(~,~) updateFields();
    center2Field.ValueChangedFcn = @(~,~) updateFilds();
    updateButton.ButtonPushedFcn = @(~,~) updateRecon();
    
%     function processGE2()
%         % Check for and read new ge2 image
%         [filePath, fileName] = selectFile(onlineDir);
%         try
%             t = 1;
%             [b, fname] = readGE2();
%         catch
%             return
%         end
% 
%         fileLabel.Text = fname; % Update file name in gui
%         t = 1;  % Increment file counter
%         
%         % Update dictionary if size of b changes
%         updateDPX();
% 
%         processImage()
%     end

    function updateDPX()
        if  ~isequal(size(b), oldbsize)
            P.N1 = size(b,1);
            P.N2 = size(b,2);
            P.mu1 = round(P.N1/2);
            P.mu2 = round(P.N2/2);
            D = dictionary2D(P);
            Df = fft2(D);
            updateDictionary();
            x_init = zeros(P.N1,P.N2,P.K);
        end
        if isempty(x_init)
            x_init = zeros(P.N1,P.N2,P.K);
        end
    end

    function updateRecon()
        try
            [b,~] = readGE2();
        catch
            return
        end

        % Update dictionary if size of b changes
        updateDPX();

        processImage();
    end
    function updateFields()
        params.lambda = lambdaField.Value;
        r1 =  rad1Field.Value;
        r2 =  rad2Field.Value;
        center(1) = center1Field.Value;
        center(2) = center2Field.Value;
    end

    % Start the timer
    function startTimer()
        start(timerObj);
        disp('Processing Started')
    end

    % Pause the timer
    function pauseTimer()
        stop(timerObj);
        disp('Processing Paused')
    end

    % Reset the data
    function resetData()
        D = [];  % Clear the dictionary
        vdf = [];  % Clear the vdf
        b = [];
        normb = 1;
        awmv_eta = [];
        awmv_2th = [];
        processedFiles = {};
        t = 0;
        cla(dictionaryAxes);
        cla(vdfAxes);
        cla(awmvPlot);
        cla(bImageAxes);
        cla(bhatImageAxes);
        fileLabel.Text = 'No file being processed';
        disp('Processing Reset')
    end

    % Timer callback function
    function timerCallback(~,~)
        if ~continueBox.Value
            pause(0.01)
            return
        end
        % Check for and read new ge2 image
        try
            [b, fname] = readH5online();
        catch
            return
        end

        fileLabel.Text = fname; % Update file name in gui
        t = t + 1;  % Increment file counter
        

        % Update dictionary if size of b changes
        updateDPX();

        processImage()
        processedFiles = [processedFiles, fname];
    end

    % processH5 callback function
    function processH5()
            
        % Read in a selected H5 image
        [fPath, fileName] = selectFile(filePath);
                
        % Check for matching output files to preload those instead
        outPattern = fullfile(outPath,'*.mat');
        outFiles = dir(outPattern);
        outFileCell = {};
        for i = 1:numel(outFiles)
            [~,subPart] = fileparts(outFiles(i).name);
            outFileCell{i} = subPart;
        end
        
        % Load solution if it exists
        if ismember(fileName, outFileCell)
            outVars = load(fullfile(outPath,[fileName,'.mat']),'X','b','P');
            x_init = outVars.X;
            b = outVars.b;
            oldbsize = size(b);
            P = outVars.P;
            D = dictionary2D(P);
            Df = fft2(D);
            t = t + 1;
            updateDictionary();
            fname = fileName;
        else
            try
                fname = fullfile(fPath,fileName);
                [b, rad, az] = loadH5polar(omInd,omInd2,fname,center,r1,r2);
                if isempty(b)
                    error('b is empty')
                end
                normb = 10/norm(b(:));
                b = b./norm(b(:))*10;

                % Need to subtract the background
                b_median = median(b(:));
                b = b-1.1*b_median;
                b(b<0) = 0;
                t = t + 1;
            catch
                disp('File not processed')
                return
            end
        end
        
        rad1Field.Value = r1;
        rad2Field.Value = r2;
        center1Field.Value = center(1);
        center2Field.Value = center(2);

        fileLabel.Text = fileName; % Update file name in gui
        
        % Update dictionary if size of b changes
        updateDPX();

        processImage()
    end

    function nextImage()
        % Read in a next H5 image
        fname = fullfile(filePath,fileName);
        
        % Check for matching output files to preload those instead
        outPattern = fullfile(outPath,'*.mat');
        outFiles = dir(outPattern);
        outFileCell = {};
        for i = 1:numel(outFiles)
            [~,subPart] = fileparts(outFiles(i).name);
            outFileCell{i} = subPart;
        end
        
        % Load solution if it exists
            if ismember(fileName, outFileCell)
                outVars = load(fullfile(outPath,[fileName,'.mat']),'X','b','P');
                x_init = outVars.X;
                b = outVars.b;
                oldbsize = size(b);
                P = outVars.P;
                D = dictionary2D(P);
                Df = fft2(D);
                updateDictionary();
                fname = fileName;
                return; 
            end
        
        try
            t = t + 1;
            [b, rad, az] = loadH5polar(omInd,omInd2,fname,center,r1,r2);
        catch
            disp('File not loaded')
            t = t - 1;
            return
        end
        % Read next ge2 image
%         try
%             t = t + 1;
%             [b, fname] = readGE2();
%         catch
%             disp('File not processed')
%             t = t - 1;
%             return
%         end

        fileLabel.Text = fileName; % Update file name in gui        
        
        % Update dictionary if size of b changes
        updateDPX();

        processImage()
    end

    % Update the dictionary display
    function updateDictionary()
        cla(dictionaryAxes);
        
        % Calculate the layout of subplots
        rows = P.K1;
        cols = P.K2;
        
        % Calculate the width and height of each dictionary entry
        entryHeight = P.N1;
        entryWidth = P.N2;
        
        % Calculate the size of the dictionary image
        dictionaryWidth = cols * entryWidth;
        dictionaryHeight = rows * entryHeight;
        
        % Create the dictionary image
        dictionaryImage = zeros(dictionaryHeight, dictionaryWidth);
        
        rowInd = 1;
        colInd = ones(P.K2,1);
        cropWidths = zeros(P.K2,1);
        fRow = 1;
        fCol = 1;

        Dshape = reshape(D,[P.N1,P.N2,P.K1,P.K2]);
        % Arrange the dictionary entries into the dictionary image
        for i = 1:P.K1
            for j = 1:P.K2
                % Get the current dictionary entry
                entry = Dshape(:,:,i,j);
                
                % Find the smallest rectangle that contains nonzero pixels
                rowSum = sum(entry,2);
                colSum = sum(entry,1);
                keepRows = rowSum>max(rowSum)/20;
                keepCols = colSum>max(colSum)/20;
                croppedImage = entry(keepRows,keepCols);
                if (i > 1)
                    cWidth = round((cropWidths(j)-size(croppedImage,2))/2);
                else
                    cropWidths(j) = size(croppedImage,2);
                    cWidth = 0;
                end

                % Copy the cropped image to the dictionary image
                dictionaryImage(rowInd:rowInd+size(croppedImage,1)-1,...
                                colInd(j)+cWidth:colInd(j)+cWidth+size(croppedImage,2)-1) = croppedImage/max(croppedImage(:));
                if (i == 1) && (j<P.K2)
                    colInd(j+1) = colInd(j) + size(croppedImage,2);
                end
                fRow = max(fRow,rowInd+size(croppedImage,1)-1);
                fCol = max(fCol,colInd(j)+size(croppedImage,2)-1);
            end
            rowInd = rowInd + size(croppedImage,1);
        end
        dictionaryImage = dictionaryImage(1:fRow,1:fCol);
        
        % Display the dictionary image using imagesc
        imagesc(dictionaryAxes, dictionaryImage);

    end

    % Update the VDF display
    function updateVDF()
        cla(vdfAxes);
        imagesc(vdfAxes, vdf);
    end

    % Update the AWMV plot
    function updateAWMVPlot()
        cla(awmvPlot);
        tArray = 1:t;
        plot(awmvPlot, tArray, awmv_eta, 'b-o', tArray, awmv_2th, 'r-x');
        legend(awmvPlot, 'AWMV_\eta', 'AWMV_{2\theta}');
    end

    % Update the b and b_hat images
    function updateImages(b_hat)
        cla(bImageAxes);
        cla(bhatImageAxes);
        
        imagesc(bImageAxes, (b));
        caxis(bImageAxes, [min(b(:)),max(b(:))]);
        title(bImageAxes,sprintf('Data %0.2f',max(b(:))/normb))

        imagesc(bhatImageAxes,(b_hat));
        caxis(bhatImageAxes, [min(b(:)),max(b(:))]);
        title(bhatImageAxes,sprintf('Reconstruction, Relative Error = %0.2f',rel_error))
    end

function [b,fname]= readH5online()
        % Check for new files in the directory

        files = dir(filePattern);
        files2 = dir(filePattern2);
        [fileNums] = cellfun(@(x)sscanf(x,'ff1_%d.h5'), {files.name});
        [sfNums,sortInd] = sort(fileNums);
        sortInd(sfNums<startFnum) = [];
        files = files(sortInd);
        files2 = files2(sortInd);
        
        % Check for matching output files to preload those instead
        outPattern = fullfile(outPath,'*.mat');
        outFiles = dir(outPattern);
        outFileCell = {};
        for i = 1:numel(outFiles)
            [~,subPart] = fileparts(outFiles(i).name);
            outFileCell{i} = subPart;
        end
        
        if isempty(files)
            if printDot
                fprintf('.')
            else
                disp('No new files')
                printDot = 1;
            end
            return;  % No new files found
        else
            printDot = 0;
                
        end
        
        % Process each file
        for i = 1:numel(files)
            fileName = files(i).name;
            
            fPath = fullfile(files(i).folder, files(i).name);
            fPath2 = fullfile(files(i).folder,files2(i).name);
            
            if ismember(fileName, processedFiles)
                continue; % Skip processing
            end
            
            % Load solution if it exists
            if ismember(fileName, outFileCell)
                outVars = load(fullfile(outPath,[fileName,'.mat']),'X','b','P');
                x_init = outVars.X;
                b = outVars.b;
                oldbsize = size(b);
                P = outVars.P;
                normb = P.normb;
                D = dictionary2D(P);
                Df = fft2(D);
                updateDictionary();
                fname = fileName;
                return; 
            end

            % Read the image file
            try
                b = loadH5polarDex(omInd,omInd2,fPath,fPath2,center,r1,r2);
                normb = 10/norm(b(:));
                b = b./norm(b(:))*10;
                
                % Need to subtract the background
                b_median = median(b(:));
                b = b-1.1*b_median;
                b(b<0) = 0;

                fname = files(i).name;
                return
            catch
                disp(['Error, skipping ',files(i).name,' for now'])
                continue
            end
        end
    end



%     function [b,fname]= readGE2online()
%         % Check for new files in the directory
%         filePattern = fullfile(onlineDir, 'ff_0*.ge2');  % Modify the pattern to match your file type
%         files = dir(filePattern);
%         
%         if isempty(files)
%             return;  % No new files found
%         end
%         
%         % Process each file
%         for i = 1:numel(files)
% 
%             if ismember(files(i).name, processedFiles)
%                 continue; % Skip processing
%             end
% 
%             filePath = files(i).folder;
%             fileName = files(i).name;
%             fPath = fullfile(files(i).folder, files(i).name);
% 
%             % Read the image file (replace with your own file reading code)
%             try
%                 b = loadGE2polar(1,fPath,center,r1,r2);
%                 b = b./norm(b(:))*10;
%                 
%                 % Need to subtract the background
%                 b_median = median(b(:));
%                 b = b-1.1*b_median;
%                 b(b<0) = 0;
% 
%                 fname = files(i).name;
%                 return
%             catch
%                 disp(['Error, skipping ',files(i).name,' for now'])
%                 continue
%             end
%         end
%     end

%     function [b,fname]= readGE2()
%         fPath = fullfile(filePath, fileName);
% 
%         % Read the image file (replace with your own file reading code)
%         try
%             [b] = loadGE2polar(t,fPath,center,r1,r2);
%             normb = 10/norm(b(:));
%             b = b./norm(b(:))*10;
%             
%             % Need to subtract the background
%             b_median = median(b(:));
%             b = b-1.1*b_median;
%             b(b<0) = 0;
% 
%             fname = fileName;
%             return
%         catch
%             disp([fileName, ' not loaded.'])
%         end
%     end

    function processImage()
        % Call the processing function (replace with your own processing code)
        [X,L] = FISTA_Circulant_gpu(Df,b,x_init,params);
        x_init = X;
%         params.L = L;
    
        Xf = fft2(X);
        b_hat = real(ifft2(Ax_gpu(Df,Xf)));
        [b_hat,X] = gather(b_hat,X);
        rel_error = norm(b(:)-b_hat(:))/norm(b(:));

        % Update the variables and GUI components
        vdf = reshape(squeeze(sum(X,[1,2])),[P.K1,P.K2]);
        [awmv_2th_t, awmv_eta_t] = computeAWMV(X,P);
        awmv_2th(t) = awmv_2th_t;
        awmv_eta(t) = awmv_eta_t;
        
        % Update GUI components
        updateVDF();
        updateAWMVPlot();
        updateImages(b_hat);
        drawnow           
        oldbsize = size(b);
        P.normb = normb;
        save(fullfile(outPath,[fileName,'.mat']),'X','b','P')
    end
    
    function [filePath, fileName] = selectFile(defaultPath)
        [fileName, filePath] = uigetfile(fullfile(defaultPath,'*.*'), 'Select a file');
        if isequal(fileName, 0) || isequal(filePath, 0)
            disp('File selection cancelled.');
            filePath = '';
            fileName = '';
        else
            disp(['Selected file: ' fullfile(filePath, fileName)]);
        end
    end
end
