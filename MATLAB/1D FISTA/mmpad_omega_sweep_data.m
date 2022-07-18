%% Alternate Omega sweep of MMPAD data
top_dir = 'E:\MMPAD_omega';
sub_dir = {'omega4','omega3','omega2'};
sweeps = [8,65;...
          15,58;...
          22,50];

for i = 0:545
    for j = 1:3
        fName = [' mmpad_img_',num2str(i),'.mat'];
        inFile = fullfile(top_dir,fName);
        load(inFile)

        polar_image = squeeze(sum(mmpad_img(sweeps(j,1):sweeps(j,2),:,:),1));

        fName = ['mmpad_img_',num2str(i+1),'.mat'];
        outFile = fullfile(top_dir,sub_dir{j},fName);
        mkdir(fullfile(top_dir,sub_dir{j}))
        save(outFile,'polar_image')
    end
end

%% Split into rings
top_dir = 'E:\MMPAD_omega';
sub_dir = {'omega4','omega3','omega2'};
ring_dir = {'ring1','ring2','ring3','ring4'};
for i = 1:546
    for j = 1:3
        for r = 1:4
            if i == 1
                mkdir(fullfile(top_dir,sub_dir{j},ring_dir{r}))
            end
            fName = ['mmpad_img_',num2str(i),'.mat'];
            inFile = fullfile(top_dir,sub_dir{j},fName);
            load(inFile)
            
            if r == 1
                polar_image = polar_image(2:262,5:40);
            elseif r == 2
                polar_image = polar_image(2:262,210:260);
            elseif r == 3
                polar_image = polar_image(2:262,268:310);
            elseif r ==4
                polar_image = polar_image(2:262,355:395);
            end
            
            outFile = fullfile(top_dir,sub_dir{j},ring_dir{r},fName);
            save(outFile,'polar_image')
            
        end
    end
end