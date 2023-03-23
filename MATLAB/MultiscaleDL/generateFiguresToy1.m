function generateFiguresToy1(topDir,outputs,suffix)
mkdir(topDir)
D = outputs.D;
X = outputs.X;
scales = outputs.scales;
N = outputs.N;
y = outputs.y;
K = outputs.K;

Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
Utotal = sum(Uarray);

AD = reSampleCustomArray(N,D,scales);
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric'));

% Show dictionary
f1 = figure;
for i = 1:Utotal
    subplot(7,7,i)
    plot(AD(:,:,i),'Linewidth',1)
 set(gca, 'XtickLabel','')
set(gca, 'FontSize', 16)
end
f1.Position = [1 100 1800 500];
saveas(f1,fullfile(topDir,['dict',suffix,'.png']))

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
saveas(f2,fullfile(topDir,['recon',suffix,'.png']))

% Recovered VDF(t)
f3 = figure;
imagesc(squeeze(sum(sum(X,1),2)))
% set(findobj(gca, 'Type', 'line'), 'LineWidth', 30)
f3.Position = [800 100 600 300];
set(gca, 'FontSize', 20)
f3.Position = [1 100 800 400];
saveas(f3,fullfile(topDir,['vdf',suffix,'.png']))

% Spatial placement of atom groups
i = 1;
for k = 1:K
    figure;
    Ui = size(scales{k},2) + i - 1;
    imagesc( squeeze(sum(X(:,:,i:Ui,:),3)) )
    i = i + Ui;
    saveas(gcf,fullfile(topDir,['X',num2str(k),suffix,'.png']))
end

%% Remove whitespace from pngs
% Get all png files in the directory
files = dir(fullfile(topDir, '*.png'));

% Loop through each file
for i = 1:length(files)
    % Get the file name and path
    fileName = files(i).name;
    filePath = fullfile(topDir, fileName);
    
    % Read the image file
    img = imread(filePath);
    
    % Convert to grayscale and invert colors
    gray = 255 * (rgb2gray(img) < 128);
    
    % Find non-zero pixels (image content)
    [row,col] = find(gray);

    % Crop the image using the bounding box
    cropped = img( min(row):max(row), min(col):max(col),:);
    
    % Save the cropped image overwriting
    imwrite(cropped, filePath);
end


end