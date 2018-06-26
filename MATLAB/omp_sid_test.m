 %% Load data 
 disp('Loading Data...')
% 311 Ring dataset
dataset = fullfile('D:','CHESS_data','al7075_311_polar_reduced');

% K = 185*5;
% X = zeros(41,2048,K);
% k = 1;
% load polar images
% for load_step = 0:4
%     fprintf('Load %i/5\n',load_step+1)
% 	for img_num = 0:184
%         str1 = sprintf('%i',load_step);
%         str2 = sprintf('%i',img_num);
%         fileName = ['polar_image_',str1,'_',str2,'.mat'];
%         fileDir = fullfile(dataset,fileName);
%         load(fileDir)
%         X(:,:,k) = polar_image;
%         k = k + 1;
%     end
% end

%% Create dictionary in SID format
disp('Creating Dictionary...')
% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;
P.sampleDims = [37,5];

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  2,    P.num_var_r).^2;

A0_stack = unshifted_basis_matrix_stack_norm2(P);
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);

m = 1;
y = [];
for i = 1:size(A0_stack,3)
    for j = 1:size(A0_stack,4)
        b = A0_stack(:,:,i,j);
        rad_width = 4*sqrt(P.var_rad(j));
        theta_width = 4*sqrt(P.var_theta(i))*(P.num_theta/(2*pi));
        b = shift2D(b,round(rad_width),round(theta_width));
        Qv(m,1) = round(2*rad_width)-1;
        Qh(m,1) = round(2*theta_width)-1;
        b_crop = b(1:Qv(m),1:Qh(m));
        m = m + 1;
        y = [y; b_crop(:)];
%         imagesc(b_crop,[0 1])
%         pause

    end
end

M = numel(Qv);
Pv = ones(M,1);
Ph = ones(M,1);

sid = initSID2D(Qv,Qh,Pv,Ph,y);
P.bSize = [16,64; 64,16; 32,32];
P.s = 400;

%% fit omp
outDir = 'D:\CHESS_data\al7075_311_norm2_omp_sid2d';

for load_step = 1
    fprintf('Load %i/5\n',load_step+1)
	for img_num = 78:184
        fprintf('Image %i/185\n',img_num+1)
        str1 = sprintf('%i',load_step);
        str2 = sprintf('%i',img_num);
        fileName = ['polar_image_',str1,'_',str2,'.mat'];
        P.load_step = load_step;
        P.img_num = img_num;
        fileDir = fullfile(dataset,fileName);
        load(fileDir)
        x_hat = saSID2Domp(polar_image, sid, P.s,'mexOMP',P.bSize,0);
        b_fit = multSID2D(sid, x_hat, 41, 2048);
        err = norm(polar_image-b_fit)/norm(polar_image);
        outName = sprintf('fista_fit_%i_%i.mat',load_step,img_num);
        save(fullfile(outDir,outName),'x_hat','polar_image','P','err')
	end
end

