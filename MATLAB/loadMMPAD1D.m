function [y,N,T] = loadMMPAD1D

% Load MMPAD subset
data_dir = 'D:\MMPAD_data_nr1\ring1_zero';
t = 1:4:100;
T = numel(t);
for i = 1:T
	load(fullfile(data_dir,['mmpad_img_',num2str(t(i)),'.mat']))
    if t == 1
        [N,M] = size(polar_image(:,134:end));
        Y = zeros(N,M,T);
        y = zeros(M,T);
    end
    Y(:,:,i) = polar_image(:,134:end)/norm(polar_image(:,134:end),'fro');
    y(:,i) = sum(Y(:,:,i),1)/norm(sum(Y(:,:,i),1),'fro');
end
[N,T] = size(y);
end

