%% Fixed Parameters
% dataset = 'D:\MMPAD_data';
% ringName = 'ring1_zero';
% ring_num  = 1;
% img_num = 134;
% prefix = 'mmpad_img';
% % Ring sampling parameters
% load(fullfile(dataset,ringName,[prefix,'_',num2str(img_num),'.mat']));
% P.set = ring_num;
% P.img = img_num;
% P.num_theta= size(polar_image,2);
% P.num_rad = size(polar_image,1);
% P.dtheta = 1;
% P.drad = 1;
% P.sampleDims = [546,1];
% 
% % Basis function variance parameters
% P.basis = 'norm2';
% P.num_var_t = 8;
% P.num_var_r = 12;
% P.var_theta = linspace(P.dtheta/2,6,P.num_var_t).^2;
% P.var_rad   = linspace(P.drad/2,  32,P.num_var_r).^2;
% 
% % Zero padding and mask
% maskRows = 129:133;
% zPad = [0,0];
% zMask = zeros(size(zeroPad(polar_image,zPad)));
% zMask(maskRows,:) = 1;
% zMask = onePad(zMask,zPad);
% [r,c] = find(zMask==1);
% zMask = [r,c];
% 
% % fista params
% params.stoppingCriterion = 1;
% params.tolerance = 1e-6;
% params.L = 1000;
% params.lambda = 0.1;
% params.gamma = 10;
% params.beta = 1.1;
% params.maxIter = 500;
% params.maxIterReg = 500;
% params.isNonnegative = 1;
% params.zeroPad = zPad;
% params.zeroMask = zMask;
% params.noBacktrack = 0;
% params.plotProgress = 0;
% P.params = params;

P.set = 1;
P.img = 134;
data_dir = 'D:\MMPAD_data\ring1_zero';
output_dir = 'D:\MMPAD_data\init_mmpad_reg_fit';
wrap_awmv_reg_FISTA_Circulant( data_dir,P,output_dir )

%% Plot results