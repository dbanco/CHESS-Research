function generateFiguresToy1zpad_center(topDir,outputs,suffix,gridDim,useMin)
if ~isempty(topDir)
    mkdir(topDir)
end
if nargin < 5
    useMin = true;
end

if useMin
    D = outputs.Dmin;
    X = outputs.Ymin;
else
    D = outputs.D;
    X = outputs.Y;
end

scales = outputs.scales;
N = outputs.N;
M = outputs.M;
y = outputs.y;
K = outputs.K;
center = (M+1)/2;
Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
Utotal = sum(Uarray);

AD = reSampleCustomArrayCenter3(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));
Yhat = unpad(Yhat,M-1,'pre');

% Show dictionary
f1 = figure;
for i = 1:Utotal
    subplot(gridDim(1),gridDim(2),i)
    plot(real(AD(:,1:M,i)),'Linewidth',2)
%     set(gca, 'XtickLabel','')
    set(gca, 'FontSize', 16)
    ylim([0,0.75])
end
f1.Position = [1 100 1800 500];
if ~isempty(topDir)
    saveas(f1,fullfile(topDir,['dict',suffix,'.png']))
end

% Recon and data regular scale
f2 = figure;
subplot(2,1,1)
imagesc(squeeze(y))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
subplot(2,1,2)
imagesc(squeeze(Yhat))
 set(gca, 'YtickLabel','')
set(gca, 'FontSize', 20)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat,'fro')/norm(y(:),'fro')))
f2.Position = [1 100 800 500];
if ~isempty(topDir)
    saveas(f2,fullfile(topDir,['recon',suffix,'.png']))
end

% Recovered VDF(t)
f3 = figure;
imagesc(squeeze(sum(sum(X,1),2)))
ylabel('\sigma')
xlabel('time')
% set(findobj(gca, 'Type', 'line'), 'LineWidth', 30)
f3.Position = [800 100 600 300];
set(gca, 'FontSize', 20)
f3.Position = [1 100 800 400];
if ~isempty(topDir)
    saveas(f3,fullfile(topDir,['vdf',suffix,'.png']))
end

% Spatial placement of atom groups
i = 1;
for k = 1:K
    figure;
    Ui = size(scales{k},2) + i - 1;
    imagesc( squeeze(sum(X(:,:,i:Ui,:),3)) )
    ylabel('shift')
    xlabel('time')
    i = i + size(scales{k},2);
    saveas(gcf,fullfile(topDir,['X',num2str(k),suffix,'.png']))
end

% Recovered X at times [1,21,31,40]
f4 = figure;
jj = 1;
T = size(X,4);
for t = [1,round(T/4),round(T/2),round(3*T/4)]
    subplot(2,2,jj)
    imagesc(squeeze(X(:,:,:,t))')
    ylabel('\sigma')
    xlabel('shift')
    title(['t=',num2str(t)])
    set(gca, 'FontSize', 18)
    colorbar()
    jj = jj + 1;
end

f4.Position = [1 100 1200 500];
if ~isempty(topDir)
    saveas(f4,fullfile(topDir,['X_times',suffix,'.png']))
end

%% Remove whitespace from pngs
% Get all png files in the directory
files = dir(fullfile(topDir, '*.png'));

% Loop through each file
for i = 1:length(files)
    % Get the file name and path
    fileName = files(i).name;
    filePath = fullfile(topDir, fileName);

    removeWhiteSpace(filePath)
end

close(f1); 
close(f2); 
close(f3);
close(f4);

end