function y = loadMMPAD1Dexp1()
topDir = 'E:\MMPAD_omega\full';
fName = 'mmpad_img_%i.mat';
skip = 4;
subset = 100;

% Load first image
f=load(fullfile(topDir,sprintf(fName,0)));
[om,eta,~] = size(f.mmpad_img);
T = 546;
sub_eta = floor(eta/2);
% Load subset of data
y = zeros(om,sub_eta,subset/skip);
j = 1;
for i = skip:skip:subset
    f=load(fullfile(topDir,sprintf(fName,i-1)));
    x = sum(f.mmpad_img,3);
    y(:,:,j) = 10*x(:,1:sub_eta)/norm(x(:,1:sub_eta));
    j = j+1;
end
end

