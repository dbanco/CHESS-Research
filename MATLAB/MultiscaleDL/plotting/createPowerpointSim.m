function createPowerpointSim(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,lambdaVals,dirStartOF,lambdaOFVals,lambdaHSVals,selected_lam_of_vec,selected_lam_s_vec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Import PowerPoint package
import mlreportgen.ppt.*;

% Create a presentation object
ppt_sim1 = Presentation(pptFile);

% Add a slide to the presentation
slide1 = add(ppt_sim1, 'Title Slide');
titleText = Paragraph(titleStr);
replace(slide1, 'Title', titleText);

NN = numel(sigmas);

for n = 2:NN
    spDir = [dirStartS,'_sig_',num2str(n)];
    hsDir = [dirStartOF,'_sig_',num2str(n)];
    j_s = find(lambdaVals == selected_lam_s_vec(n));
    j_of = find(lambdaOFVals == selected_lam_of_vec(n));
    dictPngS = sprintf("dict_j%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e.png",...
            j_s,j_of,sigmas(n),selected_lam_s_vec(n),0);
%     reconPngS = sprintf("recon_j1_sig_%0.2e_lam1_%0.2e_lam2_%0.2e.png",...
%             sigmas(n),selected_lam_s_vec(n),0);
    vdfPngS = sprintf("vdf_j%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e.png",...
            jh_s,j_of,sigmas(n),selected_lam_s_vec(n),0);
    x1PngS = sprintf("X1_j%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e.png",...
            jh_s,j_of,sigmas(n),selected_lam_s_vec(n),0);
    dictPngOF = sprintf("dict_j%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e.png",...
                j_s,j_of,sigmas(n),selected_lam_s_vec(n),selected_lam_of_vec(n));
%     reconPngOF = sprintf("recon_j%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e.png",...
%             j_of,sigmas(n),selected_lam_s_vec(n),selected_lam_of_vec(n));
    vdfPngOF = sprintf("vdf_j%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e.png",...
            j_s,j_of,sigmas(n),selected_lam_s_vec(n),selected_lam_of_vec(n));
    x1PngOF = sprintf("X1_j%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e.png",...
            j_s,j_of,sigmas(n),selected_lam_s_vec(n),selected_lam_of_vec(n));

    % Dictionaries
    slide2 = add(ppt_sim1, 'Title and Content');
    replace(slide2, 'Title', sprintf('Noise Level %i, SNR = %0.2f, Lam_of=%0.2f',...
        n,meanSNR(n),0));
    imagePath = fullfile(topDir,spDir,dictPngS);  % Specify the image file you want to insert
    img = Picture(imagePath);
    replace(slide2, 'Content', img);

    slide2 = add(ppt_sim1, 'Title and Content');
    replace(slide2, 'Title', sprintf('Noise Level %i, SNR = %0.2f, Lam_of=%0.2f',...
        n,meanSNR(n),selected_lam_of_vec(n)));
    imagePath = fullfile(topDir,hsDir,dictPngOF);  % Specify the image file you want to insert
    img = Picture(imagePath);
    replace(slide2, 'Content', img);

    % VDFs
    slide2 = add(ppt_sim1, 'Title and Content');
    replace(slide2, 'Title', sprintf('Noise Level %i, SNR = %0.2f, Lam_of=%0.2f',...
        n,meanSNR(n),0));
    imagePath = fullfile(topDir,spDir,vdfPngS);  % Specify the image file you want to insert
    img = Picture(imagePath);
    replace(slide2, 'Content', img);

    slide2 = add(ppt_sim1, 'Title and Content');
    replace(slide2, 'Title', sprintf('Noise Level %i, SNR = %0.2f, Lam_of=%0.2f',...
        n,meanSNR(n),selected_lam_of_vec(n)));
    imagePath = fullfile(topDir,hsDir,vdfPngOF);  % Specify the image file you want to insert
    img = Picture(imagePath);
    replace(slide2, 'Content', img);

    % X1
    slide2 = add(ppt_sim1, 'Title and Content');
    replace(slide2, 'Title', sprintf('Noise Level %i, SNR = %0.2f, Lam_of=%0.2f',...
        n,meanSNR(n),0));
    imagePath = fullfile(topDir,spDir,x1PngS);  % Specify the image file you want to insert
    img = Picture(imagePath);
    replace(slide2, 'Content', img);

    slide2 = add(ppt_sim1, 'Title and Content');
    replace(slide2, 'Title', sprintf('Noise Level %i, SNR = %0.2f, Lam_of=%0.2f',...
        n,meanSNR(n),selected_lam_of_vec(n)));
    imagePath = fullfile(topDir,hsDir,x1PngOF);  % Specify the image file you want to insert
    img = Picture(imagePath);
    replace(slide2, 'Content', img);

end

% Close and save the presentation
close(ppt_sim1);

end