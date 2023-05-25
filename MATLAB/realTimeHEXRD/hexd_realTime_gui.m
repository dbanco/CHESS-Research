function hexd_realTime_gui
    % Create the GUI figure
    fig = uifigure('Name', 'hexd_realTime_gui', 'Position', [10 40 1500 520]);
    
    % Create UI components
    startButton = uibutton(fig, 'Text', 'Start/Resume', 'Position', [20 20 120 40]);
    pauseButton = uibutton(fig, 'Text', 'Pause', 'Position', [160 20 120 40]);
    resetButton = uibutton(fig, 'Text', 'Reset', 'Position', [300 20 120 40]);
    
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
    title(bImageAxes,'log_{10}(Data)')
    xlabel(bImageAxes, '\eta');
    ylabel(bImageAxes, '2\theta');
    set(bImageAxes,'xtick',[],'ytick',[])

    axis(bhatImageAxes, 'tight');
    title(bhatImageAxes,'log_{10}(Reconstruction)')
    xlabel(bhatImageAxes, '\eta');
    ylabel(bhatImageAxes, '2\theta');
    set(bhatImageAxes,'xtick',[],'ytick',[])
    
    % Initialize variables
    processedFiles = {};
    onlineDir = 'C:\Users\dpqb1\Documents\Data\c103_Feb\onlineDir';
%     onlineDir = '/nfs/chess/user/dbanco/realTimeHEXRD/onlineDir';
    x_init = [];
    D = []; 
    vdf = [];
    awmv_eta = [];
    awmv_2th = [];
    t = 0;  % Number of files processed
    rel_error = 1;

    % Algoriothm setup
    center = [1025,1020];
    r1 = 430;
    r2 = 450;
    
    % Parameters
    % Basis function variance parameters
    P.K1 = 10;
    P.K2 = 15;
    P.K = P.K1*P.K2;
    P.sigma1 = linspace(0.5,  2,    P.K1);
    P.sigma2 = linspace(0.5,  20,   P.K2);
    P.basis = 'norm2';
 
    
    % fista params
    params.lambda = 0.05; % sparsity penalty
    params.L = 1000;  %
    params.beta = 2; %
    params.noBacktrack = 0;
    
    params.isNonnegative = 1; % flag to enforce nonnegativity
    
    params.stoppingCriterion = 'OBJECTIVE_VALUE';
    params.maxIter = 3;
    params.tolerance = 1e-3;
    
    params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
    params.verbose = 0;      % flag to print objective values at each iteration 
    P.params = params;

    % Set up timer
    timerObj = timer('ExecutionMode', 'fixedRate', 'Period', 1, 'TimerFcn', @timerCallback);
    
    % Button callbacks
    startButton.ButtonPushedFcn = @(~,~) startTimer();
    pauseButton.ButtonPushedFcn = @(~,~) pauseTimer();
    resetButton.ButtonPushedFcn = @(~,~) resetData();
    
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

            filePath = fullfile(files(i).folder, files(i).name);

            fileLabel.Text = files(i).name;
            % Read the image file (replace with your own file reading code)
            try
                b = loadGE2polar(1,filePath,center,r1,r2);
                b = b./norm(b(:))*10;
            catch
                disp(['Error, skipping ',files(i).name,' for now'])
                continue
            end
            
            t = t + 1;  % Increment file counter
            % b-dependent parameters
            P.N1 = size(b,1);
            P.N2 = size(b,2);
            P.mu1 = round(P.N1/2);
            P.mu2 = round(P.N2/2);

            % Construct dictionary
            D = dictionary2D(P);
            Df = fft2(D);
            updateDictionary();
            if isempty(x_init)
                x_init = zeros(P.N1,P.N2,P.K);
            end

            % Call the processing function (replace with your own processing code)
            [X,L] = FISTA_Circulant_cpu(Df,b,x_init,params);
            x_init = X;
            params.L = L;
        
            Xf = fft2(X);
            b_hat = real(ifft2(Ax_cpu(Df,Xf)));
            rel_error = norm(b(:)-b_hat(:))/norm(b(:));

            % Update the variables and GUI components
            vdf = reshape(sum(X,[1,2]),[P.K1,P.K2]);
            [awmv_2th_t, awmv_eta_t] = computeAWMV(X,P);
            awmv_2th(t) = awmv_2th_t;
            awmv_eta(t) = awmv_eta_t;
            processedFiles = [processedFiles, files(i).name]
            
            % Update GUI components
            updateVDF();
            updateAWMVPlot();
            updateImages(b, b_hat);
            drawnow
        end
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
%                     cWidth = 0;
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
    function updateImages(b, b_hat)
        cla(bImageAxes);
        cla(bhatImageAxes);
        
        imagesc(bImageAxes, log10(b));
        clim(bImageAxes, [min(log10(b(:))),max(log10(b(:)))]);

        imagesc(bhatImageAxes, log10(b_hat));
        clim(bhatImageAxes, [min(log10(b(:))),max(log10(b(:)))]);
        title(bhatImageAxes,sprintf('log_{10}(Reconstruction), Relative Error = %0.2f',rel_error))
    end
end
