function y = loadMMPAD1D(ring,interp)
% Load MMPAD subset
topDir = ['D:\MMPAD_data_nr1\ring',num2str(ring),'_zero'];
fName = 'mmpad_img_%i.mat';
% Load first image
f=load(fullfile(topDir,sprintf(fName,1)));
[theta,eta] = size(f.polar_image);

T = 546;
sub_theta = 1:theta;
sub_eta = 1:eta;
sub_T = 200;
t_step = 1;

% Load subset of data
y = zeros(numel(sub_eta),sub_T);
j = 1;
for i = t_step:t_step:sub_T
    f=load(fullfile(topDir,sprintf(fName,i)));
    
    % Sum over 2theta
    x = sum(f.polar_image(sub_theta,sub_eta),1);

    % Interpolate missing region
    if interp
        x(130:134) = (1:5)/6*(x(135)-x(129)) + x(129);
    end
    y(:,j) = x/norm(x);

    j = j+1;
end
end

