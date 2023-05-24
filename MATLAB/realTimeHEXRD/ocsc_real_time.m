clear
dbstop if error
addpath('basic_tool'); 
addpath('OCSC');
addpath('mtimesx');%**
%% set para
K = [5];
psf_s=11;                                                                                                      
psf_radius = floor( psf_s/2 );
precS = 1;
use_gpu = 0;

%% load data
Trange = 1:10;
center = [1025,1020];
r1 = 430;
r2 = 450;
factor = 1e-3;
onlineDir = 'C:\Users\dpqb1\Documents\Data\c103-90-ff-1\1\ff\onlineDir';
fname = 'polar_image_1.mat';
load(fullfile(onlineDir,fname));
padB = padarray(b, [psf_radius, psf_radius, 0], 0, 'both');
PARA= auto_para(K,psf_s,b,'all',1e-3,precS,use_gpu);
%% run
gotop