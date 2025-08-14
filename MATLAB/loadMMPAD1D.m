function y = loadMMPAD1D(ring,interp,topDir)
% Load MMPAD subset
if nargin < 3
    topDir = 'E:\Data\mmpad';
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
    if strcmp(interp, 'linear')
        x(129:134) = (1:6)/7*(x(135)-x(128)) + x(128);
    end
   
    y(:,j) = x;
    j = j+1;
end

if strcmp(interp, 'copy-shift')
        y(129,1:115) = y(134,6:120);
        y(130,1:116) = y(134,5:120);
        y(131,1:117) = y(134,4:120);  
        y(132,1:118) = y(134,3:120);
        y(133,1:119) = y(134,2:120);
  
        y(129,116:120) = y(129,115);
        y(130,117:120) = y(130,116);
        y(131,118:120) = y(131,117);  
        y(132,119:120) = y(132,118);
        y(133,120) = y(133,119);
end

for j = 1:size(y,2)
    y(:,j) = y(:,j)/norm(y(:,j));
end