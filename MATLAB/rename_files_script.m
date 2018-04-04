% Load file and save with new name

direc = 'E:\CHESS_data\al7075_311_polar_fit2\';
direcNew = 'E:\CHESS_data\al7075_311_polar_fit3\';
baseName = 'fista_fit_%i_%i_task_%i.mat';
baseNameNew = 'fista_fit_%i_%i.mat';

task_num = 0;
for load_step = 0:4
    for img_num = 0:204
        try
            fName = sprintf(baseName,load_step,img_num,task_num);
            load(fullfile(direc,fName))
            P.sampleDims = [41,5];
            
            newName = sprintf(baseNameNew,load_step,img_num);
            save(fullfile(direcNew,newName),'err','P','polar_image','x_hat')
        catch
            fName = sprintf(baseNameNew,load_step,img_num);
            load(fullfile(direc,fName))
            P.sampleDims = [41,5];
            newName = sprintf(baseNameNew,load_step,img_num);
            save(fullfile(direcNew,newName),'err','P','polar_image','x_hat')
        end
        task_num = task_num + 1
    end
end