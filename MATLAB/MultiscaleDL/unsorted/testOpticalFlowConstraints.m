% Load outputs
outputDir = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_200exact_sig_1';

fname = cell(2,1);
% sArray = cell(3,1);
% uArray = cell(3,1);
% vArray = cell(3,1);
% of_conArray = cell(3,1);

fName{1} = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_200exact_sig_1\output_j1_sig_0.00e+00_lam1_0.00e+00_lam2_0.00e+00.mat';
fName{2} = 'C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_initD0_200_sig_1\output_j1_sig_0.00e+00_lam1_1.50e-01_lam2_0.00e+00.mat';

K = 2;
smoothness = [1e-8 1 10];
maxIt = 1;

fprintf(    'Dataset   smoothness     OF1     OF2     OF3     OF4    \n')
outString = '%i        %2.2f          %2.2f   %2.2f   %2.2f   %2.2f  \n';

close all
  
for i = 1:3
    for dataNum = 1:2
        load(fName{dataNum})
        data = squeeze(outputs.X);
        data = convn(data,ones(3,3),'same');
        [N, KJ, T] = size(data);
        
        % Blob data first
        blobData = convn(data,ones(3,3),'same');

        [u,v,Fy,Fx,Ft]  = computeHornSchunkDict(blobData,K,smoothness(i),maxIt);
%         uArray{exNum} = u;
%         vArray{exNum} = v;
        Jterm = [];
        of_con = opticalFlowOp(data,u,v,1,0,Jterm);
%         out = norm(vec(of_con));
        out = sum(abs(vec(of_con(:))));

        Jterm = Ft == 0;
        of_con2 = opticalFlowOp(data,u,v,1,0,Jterm);
%         out2 = norm(vec(of_con2));
        out2 = sum(abs(vec(of_con2(:))));

        Jterm = [];
        of_con3 = opticalFlowOp(blobData,u,v,1,0,Jterm);
%         out3 = norm(vec(of_con3));
        out3 = sum(abs(vec(of_con3(:))));
        
        Jterm = Ft == 0;
        of_con4 = opticalFlowOp(blobData,u,v,1,0,Jterm);
%         out4 = norm(vec(of_con4));
        out4 = sum(abs(vec(of_con4(:))));
    
        fprintf(outString,dataNum,smoothness(i),out,out2,out3,out4)
    end
end
