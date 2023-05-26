function hexd_realTime_gui
    % Create the GUI figure
    fig = uifigure('Name', 'hexd_realTime_gui', 'Position', [10 40 1500 520]);
    
    % Create UI components
    startButton = uibutton(fig, 'Text', 'Start/Resume', 'Position', [20 20 120 40]);
    pauseButton = uibutton(fig, 'Text', 'Pause', 'Position', [160 20 120 40]);
    resetButton = uibutton(fig, 'Text', 'Reset', 'Position', [300 20 120 40]);
    processH5Button = uibutton(fig, 'Text', 'Read H5', 'Position', [500 20 120 40]);
    nextImageButton = uibutton(fig, 'Text', 'Next Image', 'Position', [640 20 120 40]);
    
    fileLabel = uilabel(fig, 'Text', 'No file being processed', 'Position', [20 495 400 30]);
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
    
    % Initialize variables
    filePath = '/nfs/chess/raw/2022-3/id1a3/miller-3528-a/ff-c103-90-s2-1/2/ff';
    fileName = [];
    processedFiles = {};

    onlineDir = '/nfs/chess/raw/2023-1/id1a3/miller-3528-b/c103-90-s2-4/3/ff';
%     onlineDir = 'C:\Users\dpqb1\Documents\Data\c103_Feb\onlineDir';

    x_init = [];
    D = []; 
    Df = [];
    vdf = [];
    awmv_eta = [];
    awmv_2th = [];
    t = 0;  % Number of files processed
    rel_error = 1;
    oldbsize = [];

    % Algorithm setup for ge2
    center = [1025,1020];
    r1 = 425;
    r2 = 445;
    
    % Parameters
    % Basis function variance parameters
    P.K1 = 10;
    P.K2 = 15;
    P.K = P.K1*P.K2;
    P.sigma1 = linspace(0.5,  4,    P.K1);
    P.sigma2 = linspace(0.5,  12,   P.K2);
    P.basis = 'norm2';
 
    
    % fista params
    params.lambda = 1e-3; % sparsity penalty
    params.L = 1000;  %
    params.beta = 2; %
    params.noBacktrack = 0;
    
    params.isNonnegative = 1; % flag to enforce nonnegativity
    
    params.stoppingCriterion = 'OBJECTIVE_VALUE';
    params.maxIter = 100;
    params.tolerance = 1e-3;
    
    params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
    params.verbose = 1;      % flag to print objective values at each iteration 
    P.params = params;

    % Set up timer
    timerObj = timer('ExecutionMode', 'fixedRate', 'Period', 1, 'TimerFcn', @timerCallback);
    
    % Button callbacks
    startButton.ButtonPushedFcn = @(~,~) startTimer();
    pauseButton.ButtonPushedFcn = @(~,~) pauseTimer();
    resetButton.ButtonPushedFcn = @(~,~) resetData();
    processH5Button.ButtonPushedFcn = @(~,~) processH5();
    nextImageButton.ButtonPushedFcn = @(~,~) nextImage();
    
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
        awmv_eta = [];
        awmv_2th = [];
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
            
        % Check for and read new ge2 image
        try
            [b,rad,az, fname] = readGE2();
        catch
            return
        end

        fileLabel.Text = fname; % Update file name in gui
        t = t + 1;  % Increment file counter
        
        % Update dictionary if size of b changes
        if  ~isequal(size(b), oldbsize)
            P.N1 = size(b,1);
            P.N2 = size(b,2);
            P.mu1 = round(P.N1/2);
            P.mu2 = round(P.N2/2);
            D = dictionary2D(P);
            Df = fft2(D);
            updateDictionary();
        end

        if isempty(x_init)
            x_init = zeros(P.N1,P.N2,P.K);
        end

        processImage(b,rad,az)
        processedFiles = [processedFiles, fname];
    end

    % processH5 callback function
    function processH5()
            
        % Read in a selected H5 image
        [filePath, fileName] = selectFile();
        try
            img = squeeze(h5read( fullfile(filePath,fileName),...
                                '/imageseries/images', ...
                                [1 1 1],[3072 3888 1] ));
        catch
            disp('File not processed')
            return
        end

        fileLabel.Text = fileName; % Update file name in gui
        t = t + 1;  % Increment file counter
        
        % Update dictionary if size of b changes
        if  ~isequal(size(b), oldbsize)
            P.N1 = size(b,1);
            P.N2 = size(b,2);
            P.mu1 = round(P.N1/2);
            P.mu2 = round(P.N2/2);
            D = dictionary2D(P);
            Df = fft2(D);
            updateDictionary();
        end

        if isempty(x_init)
            x_init = zeros(P.N1,P.N2,P.K);
        end

        processImage(b,rad,az)
    end

    function nextImage()
        % Read in a next H5 image
        try
            t = t + 1;  % Increment file counter
            img = squeeze(h5read( fullfile(filePath,fileName),...
                                '/imageseries/images', ...
                                [1 1 t],[3072 3888 t] ));
            % Still need to extract a ring to polar coordinates
        catch
            disp('File not processed')
            t = t - 1;
            return
        end

        fileLabel.Text = fileName; % Update file name in gui
        
        % Update dictionary if size of b changes
        if  ~isequal(size(b), oldbsize)
            P.N1 = size(b,1);
            P.N2 = size(b,2);
            P.mu1 = round(P.N1/2);
            P.mu2 = round(P.N2/2);
            D = dictionary2D(P);
            Df = fft2(D);
            updateDictionary();
        end

        if isempty(x_init)
            x_init = zeros(P.N1,P.N2,P.K);
        end

        processImage(b,rad,az)
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
    function updateImages(b, b_hat,rad,az)
        cla(bImageAxes);
        cla(bhatImageAxes);
        
        imagesc(bImageAxes,az,rad, (b));
        caxis(bImageAxes, [min(b(:)),max(b(:))]);

        imagesc(bhatImageAxes, az,rad,(b_hat));
        caxis(bhatImageAxes, [min(b(:)),max(b(:))]);
        title(bhatImageAxes,sprintf('Reconstruction, Relative Error = %0.2f',rel_error))
    end

    function [b,rad,az,fname]= readGE2()
        % Check for new files in the directory
        filePattern = fullfile(onlineDir, 'ff_0*.ge2');  % Modify the pattern to match your file type
        files = dir(filePattern);
        
        if isempty(files)
            return;  % No new files found
        end
        
        % Process each file
        for i = 1:numel(files)

            if ismember(files(i).name, processedFiles)
                continue; % Skip processing
            end

            fPath = fullfile(files(i).folder, files(i).name);

            % Read the image file (replace with your own file reading code)
            try
                [b,rad,az] = loadGE2polar(1,fPath,center,r1,r2);
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

    function processImage(b,rad,az)
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
        updateImages(b, b_hat,rad,az);
        drawnow           
        oldbsize = size(b);
    end
    
    function [filePath, fileName] = selectFile()
        [fileName, filePath] = uigetfile('*.*', 'Select a file');
        if isequal(fileName, 0) || isequal(filePath, 0)
            disp('File selection cancelled.');
            filePath = '';
            fileName = '';
        else
            disp(['Selected file: ' fullfile(filePath, fileName)]);
        end
    end
end
