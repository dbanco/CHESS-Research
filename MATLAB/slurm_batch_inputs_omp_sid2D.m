addpath(genpath('./'))
% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar_reduced');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','al7075_311_norm2_omp_sid2D_1');

mkdir(outputdir)
% Function
funcName = 'wrap_norm2_omp_sid2D';

jobDir = fullfile('/cluster','home','dbanco02','job_al7075_311_norm2_omp_sid2D_1');
mkdir(jobDir)

%% Fixed Parameters

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

% Prepate SID2D
A0_stack = unshifted_basis_matrix_stack_norm2(P);

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

% omp_sid2D params
P.sid = sid;
P.bSize = [16,64; 64,16; 32,32];
P.s = 200;

%% Parameters to vary
img_nums = 0:184;
load_steps = 0:4;

k = 0;
for load_step = load_steps
    for img = img_nums
        P.img = img;
        P.load_step = load_step;

        varin = {dataset,P,outputdir};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

slurm_write_bash(k-1,jobDir,'full_batch_script.sh','0-924')
