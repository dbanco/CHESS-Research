function realTimeCDL()
    % Create a figure and axes
    fig = figure('Position', [200, 200, 400, 675]);
    dictAxis = axes('Parent', fig, 'Position', [0.1, 0.75, 0.8, 0.25]);
    solutionAxis = axes('Parent', fig, 'Position', [0.1, 0.45, 0.8, 0.25]);
    dataAxis = axes('Parent', fig, 'Position', [0.1, 0.15, 0.8, 0.25]);
    t1 = title(solutionAxis, 'Recon:');
    t2 = title(dataAxis, 'Data:');
    set(t1,'Interpreter','none')
    set(t2,'Interpreter','none')
    % Create "Run" button
    runBtn = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
        'String', 'Start OCDL', ...
        'Position', [150, 10, 100, 30], ...
        'Callback', @runCallback);

    % Create "Stop" button
    stopBtn = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
        'String', 'Stop', ...
        'Position', [260, 10, 100, 30], ...
        'Callback', @stopCallback);

    % Create horizontal scrollbar for the solution axis
    scrollbar = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Position', [80, 60, 240, 20], ...
        'Callback', @scrollCallback);

    % Timer object to control the function execution
    t = timer('ExecutionMode', 'fixedRate', 'Period', 0.1, 'TimerFcn', @timerCallback);
    
    % Flag to indicate if the function is running
    isRunning = false;

    % Initialize parameters
    onlineDir = 'C:\Users\dpqb1\Documents\Data\c103_Feb\onlineDir';
    dbstop if error
    addpath('basic_tool'); 
    addpath('OCSC');
    addpath('mtimesx');%**
    
    % set para
    K = [6];
    psf_s=11;                                                                                                      
    psf_radius = floor( psf_s/2 );
    precS = 1;
    use_gpu = 0;
    b = zeros(21,2827);
    para= auto_para(K,psf_s,b,'all',1e-3,precS,use_gpu);

    A_h = [];
    B_h=[];
    d_small = init_dic(para);
    %d_small = randn(para.size_k);
    if (para.precS ==1)
        d_small = single(d_small);
    end
    if (para.gpu ==1)
        d_small = gpuArray(d_small);
    end   
    d = dsmall2d(d_small,para);
    d_hat = fft2(d); 

    s=[];
    y=[];
    scale = 1;
    s_i = 0;
    processedFiles = {};

    % Callback function for "Run" button
    function timerCallback(~, ~)
        % Reset stop flag
        if isRunning
            % Check online Dir for new .mat files
            listing = dir([onlineDir,'/polar_image*.mat']);
            process = 0;
            for f_i = 1:numel(listing)
                process = 1;
                for f_ii = 1:numel(processedFiles)
                    if strcmp(listing(f_i).name,processedFiles{f_ii})
                        process = 0;
                        break
                    end
                end
                % After checking the full list process the file
                if process
                    break
                end
            end
            % If new load and process
            if process
                % Load data file
                load(fullfile(listing(f_i).folder,listing(f_i).name),'b')
                Mtb = padarray(b, [para.psf_radius, para.psf_radius, 0], 0, 'both');
%                 mul_heur = 50;
%                 gamma_heuristic = mul_heur* 1/max(b(:));
% %                 max(b(:))
%                 para.rho_D = 1;%gamma_heuristic;
%                 para.rho_Z = 1;%gamma_heuristic;
                if (para.precS ==1)
                    Mtb = single(Mtb);
                end
                if (para.gpu ==1)
                    Mtb = gpuArray(Mtb);
                end  
                b_hat = fft2(Mtb);
                
                % Add new file to list of processed files
                s_i = s_i + 1;
                processedFiles{s_i} = listing(f_i).name;
                
                disp(['Processing: ',listing(f_i).name])
                % Expand coefficient size

                temp_b = b;
                temp_b_hat = b_hat;
                %% 1.pre-process Z
                t_Z = tic;%~~~~!!!
                [stat_Z] = precompute_H_hat_Z(d_hat, para);    
                %% 2.update Z
                [z_si,z_hat_si] = updateZ_ocsc(temp_b_hat,para,d_hat,stat_Z);
                timeZ = toc(t_Z);
                objZ = objective_online(z_hat_si,d_hat, temp_b_hat,para );
               if strcmp( para.verbose, 'all')
                   if (mod(s_i,scale)==0)
                        [ps] = eval_psnr(d_hat, z_hat_si,temp_b,para); 
                        fprintf('Z: no.img: %d, obj: %2.2f, psnr: %2.2f\n', s_i,objZ,ps)
                    end 
               end
               save(fullfile(listing(f_i).folder,'outputs',['coefs_',listing(f_i).name]),'z_hat_si')
                clear stat_Z           
                %% 1.pre-process D
                if isempty(A_h)
                    init_or_not = 1;
                else
                    init_or_not = 2;
                end
                t_D =tic;
                if para.gpu ==1
                    [A_h,B_h] = hist_ocsc_gpu(temp_b_hat,z_hat_si, para, A_h,B_h,init_or_not);
                else
                    [A_h,B_h] = hist_ocsc_cpu(temp_b_hat,z_hat_si, para, A_h,B_h,init_or_not);
                end
                %% 2.update D
                [d,d_hat,s,y] = updateD_ocsc(para,A_h,B_h,s,y,d_hat);    
                timeD =toc(t_D);
                d_curr = d2dsmall(d,para);
%                 if strcmp( para.verbose, 'all')
%                     if para.gpu==1
%                         d_show = gather(d_curr);
%                         show_dic(d_show,para,0,0);
%                     else
%                         show_dic(d_curr,para,0,0); 
%                     end
%                 end
                y_hat = real(ifft2(sum(d_hat.*z_hat_si,3)));
                y_hat = y_hat(1 + para.psf_radius:end - para.psf_radius,...
                              1 + para.psf_radius:end - para.psf_radius,:);

                idx1 = 1;
                idx2 = 50;
                maxIdx = size(y_hat,2);
                
                % Update the solution axis with the latest solution
                axes(solutionAxis);
                cla(solutionAxis);
                imagesc(y_hat)
                t1 = title(solutionAxis, ['Recon: ',listing(f_i).name]);
                set(t1,'Interpreter','none')

                axes(dataAxis);
                cla(dataAxis);
                imagesc(b)
                t2 = title(dataAxis, ['Data: ',listing(f_i).name]);
                set(t2,'Interpreter','none')
                
%                 [PSNR,RMSE] = my_psnr(b,y_hat)

                % Update the scrollbar limits
                scrollbar.Min = 1;
                scrollbar.Max = maxIdx-50;
                scrollbar.Value = 1;

                % Update the dict axis with the dictionary being learned
                axes(dictAxis);
                cla(dictAxis);
                show_dic(d_curr,para,0,0); 
                title(dictAxis, ['Dictionary: ',listing(f_i).name]);

            end
        end
    end
    function runCallback(~, ~)
        if ~isRunning
            % Start the timer and set the flag
            start(t);
            isRunning = true;
            disp('Processing started')
        end
    end

    % Callback function for "Stop" button
    function stopCallback(~, ~)
        % Set stop flag to true
        if isRunning
            stop(t);
            isRunning = false;
    
            % Save and plot outputs
            out.d = d;
            out.d_hat = d_hat;
            out.A_h = A_h;
            out.B_h = B_h;
            
            save(fullfile(onlineDir,'outputs','dict_out.mat'),'out')
            disp('Final Dictionary Saved')
        end
    end

    % Callback function for the scrollbar
    function scrollCallback(~, event)
        val = event.Source.Value;
        xlim(solutionAxis, [1+val, 50 + val]);
        xlim(dataAxis, [1+val, 50 + val]);
    end
end

%                     show_dic(d_curr,para,0,0); 
%                 end
%             end
%         end
%     end
%     disp('Stopped')
%     end
% end
% 
