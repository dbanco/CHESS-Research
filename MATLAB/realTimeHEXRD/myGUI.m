function myGUI()
    % Create a figure and axes
    fig = figure('Position', [200, 200, 400, 400]);
    dictAxis = axes('Parent', fig, 'Position', [0.1, 0.55, 0.8, 0.4]);
    solutionAxis = axes('Parent', fig, 'Position', [0.1, 0.1, 0.8, 0.4]);

    % Create "Run" button
    runBtn = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
        'String', 'Run Function', ...
        'Position', [150, 10, 100, 30], ...
        'Callback', @runCallback);

    % Create "Stop" button
    stopBtn = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
        'String', 'Stop Function', ...
        'Position', [260, 10, 100, 30], ...
        'Callback', @stopCallback);

    % Flag to control the while loop
    stopFlag = false;

    % Callback function for "Run" button
    function runCallback(~, ~)
        % Reset stop flag
        stopFlag = false;

        % Call your function here
        yourFunction(stopFlag);
    end

    % Callback function for "Stop" button
    function stopCallback(~, ~)
        % Set stop flag to true
        stopFlag = true;
    end

    % Your function that will be executed when the "Run" button is pressed
    function yourFunction(stopFlag)
        % Your code here
        % Replace this with your actual function implementation

        x = linspace(0, 2*pi, 100);
        y = sin(x);
        result = [x; y];

        % Update the solution axis with the latest solution
        cla(solutionAxis);
        plot(solutionAxis, result);
        xlabel(solutionAxis, 'X');
        ylabel(solutionAxis, 'Y');
        title(solutionAxis, 'Latest Solution');

        % Update the dict axis with the dictionary being learned
        cla(dictAxis);
        imagesc(rand(40,40))
        xlabel(dictAxis, 'X');
        ylabel(dictAxis, 'Y');
        title(dictAxis, 'Dictionary');

        % Example of a while loop that can be interrupted by the "Stop" button
        idx = 1;
        while ~stopFlag && idx <= length(x)
            % Perform computations or any time-consuming task here

            % Update the solution axis with the current progress
            cla(solutionAxis);
            plot(solutionAxis, x(1:idx), y(1:idx));
            xlabel(solutionAxis, 'X');
            ylabel(solutionAxis, 'Y');
            title(solutionAxis, 'Latest Solution');

            drawnow;  % Update the figure window

            idx = idx + 1;
            pause(0.1);  % Just for demonstration purposes, replace with your actual computation time
        end
    end
end
