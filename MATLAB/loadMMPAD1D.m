function y = loadMMPAD1D(ring,interp,topDir)
% Load MMPAD subset
if nargin < 3
    topDir = ['E:\Data\mmpad\ring',num2str(ring),'_zero'];
end

fName = 'mmpad_img_%i.mat';
% Load first image
f=load(fullfile(topDir,sprintf('ring%i_zero',ring),sprintf(fName,1)));
[theta,eta] = size(f.polar_image);

T = 546;
sub_theta = 1:theta;
sub_eta = 1:eta;
sub_T = 120;
t_step = 1;

% Load subset of data
y = zeros(numel(sub_eta),sub_T);
j = 1;
for i = t_step:t_step:sub_T
    f=load(fullfile(topDir,sprintf('ring%i_zero',ring),sprintf(fName,i)));
    
    % Sum over 2theta
    x = sum(f.polar_image(sub_theta,sub_eta),1);

    % Interpolate missing region
    if interp
        x(129:134) = (1:6)/7*(x(135)-x(128)) + x(128);
    end
    y(:,j) = x/norm(x);

    j = j+1;
end
end

