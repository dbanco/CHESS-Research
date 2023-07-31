function postAnalysisGUI
    figure('Name', 'Post Analysis GUI', 'Position', [100 100 1500 600]);

    outDir = '/nfs/chess/user/dbanco/outputs/ti-2-tension_ring1/';

    % Get a list of .mat files in the directory
    matFiles = dir(fullfile(outDir, '*.mat'));
    N = numel(matFiles);
    K = 3;
    
    % Initialize variables
    currentFileIdx = 1;
    bufferIdx = 1;
    
    buffer(1).b
    b_hat = cell(K,1);
    b = cell(K,1);
    X = cell(K,1);
    P = cell(K,1);
    vdf = cell(K,1);
    rel_err = zeros(K,1);
    awmv_eta = zeros(1, K);
    awmv_2th = zeros(1, K);
    
    for i = currentFileIdx:currentFileIdx+K-1
        fprintf('Loading %i\n',i)
        data = load(fullfile(outDir, matFiles(i).name));
        % Update the variables and GUI components
        if i == 1
            D = dictionary2D(data.P);
            Df = fft2(D);
        end
        buffer(i).X = data.X;
        buffer(i).b = data.b;
        Xf = fft2(data.X);
        buffer(i).b_hat = real(ifft2(Ax_cpu(Df,Xf)));
        buffer.P(i) = data.P;
        buffer.vdf(i) = reshape(squeeze(sum(X{i},[1,2])),[P{i}.K1,P{i}.K2]);
        [awmv_2th_t, awmv_eta_t] = computeAWMV(X{i},P{i});
        awmv_2th(i) = awmv_2th_t;
        awmv_eta(i) = awmv_eta_t;
        rel_err(i) = norm(b{i}(:)-b_hat{i}(:))/norm(b{i}(:));
    end
    
    % Create axes for b_hat
    ax1 = axes('Position', [0.05 0.05 0.9 0.2]);
    ax2 = axes('Position', [0.05 0.4 0.9 0.2]);
    ax3 = axes('Position', [0.7 0.7 0.22 0.22]);
    ax4 = axes('Position', [0.4 0.7 0.25 0.25]);
    
    updatePlots();
   
    
    % Define keyboard callback function
    set(gcf, 'KeyPressFcn', @keyboardCallback);

    % Keyboard callback function
    function keyboardCallback(~, event)
        if strcmp(event.Key, 'x')
            if currentFileIdx < numel(matFiles)
                currentFileIdx = currentFileIdx + 1;
                bufferIdx = bufferIdx + 1;
                if bufferIdx == 3
                    bufferIdx = 2;
                    % Clear 1st slot
                    buffer(1) = [];
                    % Load next into 3rd buffer spot
                    loadData(currentFileIdx+1)
                end
                updatePlots();
            end
                
        elseif strcmp(event.Key, 'z')
            if currentFileIdx > 1
                currentFileIdx = currentFileIdx - 1;
                bufferIdx = bufferIdx - 1;
                if bufferIdx == 1
                    
                    if currentFileIdx > 1
                        bufferIdx = 2;
                        % Clear 3rd slot

                        % Load next lower into 1st slot if currentFileIdx > 1
                        buffer(3) = buffer(2);
                        buffer(2) = buffer(1);
                    end
                    
                updatePlots();
                end
            end
        end
    end

    function updatePlots()
        
        % b_hat
        imagesc(ax1,buffer(bufferIdx).b_hat);
        title(ax1,sprintf('b_{hat}, Rel Error = %0.2f',rel_err(currentFileIdx)));
        
        % b
        imagesc(ax2,buffer(bufferIdx).b);
        title(ax2,['b: ',matFiles(currentFileIdx).name])
        
        % awmv
        axes(ax3)
        cla(ax3)
        hold off
        plot(awmv_eta, 'bo-');
        plot(awmv_2th, 'rx-');
        title('AWMV');
        legend('AWMV_\eta', 'AWMV_{2\theta}');
        hold on;
        plot(currentFileIdx, awmv_eta(currentFileIdx), 'gs', 'MarkerSize', 10);
        plot(currentFileIdx, awmv_2th(currentFileIdx), 'gs', 'MarkerSize', 10);
        hold off;
        
        % vdf
        axes(ax4)
        imagesc(ax4,buffer(bufferIdx).vdf)
        title('VDF')
    end

    function loadData(idx)
        data = load(fullfile(outDir, matFiles(idx).name));
        buffer(idx).X = data.X;
        buffer(idx).b = data.b;
        Xf = fft2(data.X);
        buffer(idx).b_hat = real(ifft2(Ax_cpu(Df,Xf)));
        buffer.P(idx) = data.P;
        buffer.vdf(i) = reshape(squeeze(sum(X{i},[1,2])),[P{i}.K1,P{i}.K2]);
        [awmv_2th_t, awmv_eta_t] = computeAWMV(X{i},P{i});
        awmv_2th(i) = awmv_2th_t;
        awmv_eta(i) = awmv_eta_t;
        rel_err(i) = norm(b{i}(:)-b_hat{i}(:))/norm(b{i}(:));
end
