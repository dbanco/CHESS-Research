%% Processes data to be loaded into spread_visualizer
spreadDir = fullfile('D:\CHESS_data\spread_results');
outFile = 'spread_311_norm2_omp_sid2d_1.mat';
result_path = fullfile('D:\CHESS_data','al7075_311_norm2_omp_sid2d');
fileInfo = dir(fullfile(result_path,'fista_fit*.mat'));

%% Create SID dictionary variance values 
% disp('Creating Dictionary...')
% 
% % Ring sampling parameters
% P.ring_width = 20;
% P.num_theta= 2048;
% P.num_rad = 2*P.ring_width+1;
% P.dtheta = 2*pi/P.num_theta;
% P.drad = 1;
% P.sampleDims = [37,5];
% 
% % Basis function variance parameters
% P.num_var_t = 15;
% P.num_var_r = 10;
% P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
% P.var_rad   = linspace(P.drad,  2,    P.num_var_r).^2;
% 
% m = 1;
% y_var_theta = [];
% y_var_rad = [];
% for i = 1:P.num_var_t
%     for j = 1:P.num_var_r
% 
%         rad_width = 4*sqrt(P.var_rad(j));
%         theta_width = 4*sqrt(P.var_theta(i))*(P.num_theta/(2*pi));
% 
%         Qv(m,1) = round(2*rad_width)-1;
%         Qh(m,1) = round(2*theta_width)-1;
%         b_var_rad = P.var_rad(j)*ones(Qv(m),Qh(m));
%         b_var_theta = P.var_theta(i)*ones(Qv(m),Qh(m));
%         y_var_rad = [y_var_rad; b_var_rad(:)];
%         y_var_theta = [y_var_theta; b_var_theta(:)];
%     end
% end
% 
% M = numel(Qv);
% Pv = ones(M,1);
% Ph = ones(M,1);

%% Do spread analysis

for i = 1:numel(fileInfo)
    load(fullfile(fileInfo(i).folder,fileInfo(i).name))
    if(i == 1)
        var_signal = zeros(P.num_var_t,P.num_var_r,5,37,5);
        rel_error = zeros(5,37,5);
        sparsity = zeros(5,37,5);
    end
    disp(['Load: ' num2str(P.load_step) '   Image: ' num2str(P.img)])
    idx1 = floor(P.img/5)+1;
    idx2 = mod(P.img,5)+1;
    sparsity(P.load_step+1,idx1,idx2) = sum(x_hat(:)>0);
    rel_error(P.load_step+1,idx1,idx2) = err(end);
    N_last = 1;
    for k = 1:P.num_var_t
        for j = 1:P.num_var_r
            N_coefs = (P.num_theta-Qh(k)+1)*(P.num_rad-Qv(j)+1);
            var_signal(k,j,P.load_step+1,idx1,idx2) = sum(x_hat(N_last:(N_coefs + N_last-1)));
            N_last = N_last + N_coefs;
        end
    end
end

sparsity = flip(sparsity,2);
rel_error = flip(rel_error,2);
var_signal = flip(var_signal,4);

save(fullfile(spreadDir,sprintf(outFile,P.params.lambda)),'var_signal','rel_error','sparsity','P')
        