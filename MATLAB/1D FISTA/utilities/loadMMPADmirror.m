function B = loadMMPADmirror(dataDir,Time,P)
%UNTITLED Summary of this function goes here
% Load data
B = zeros(P.num_theta,numel(Time));
for t = Time
    load(fullfile(dataDir,[P.prefix,'_',num2str(t),'.mat']))
    [n1,n2] = size(polar_image);
    if n2>n1
        polar_image = polar_image';
    end
    b = P.dataScale*sum(polar_image,2);
    % Mirror data
    nn = numel(b);
    pad1 = floor(nn/2);
    pad2 = ceil(nn/2);
    N = nn + pad1 + pad2;
    b_mirror = zeros(N,1);
    b_mirror((pad1+1):(pad1+nn)) = b;
    b_mirror((1+N-pad2):N) = flip(b((nn-pad2+1):nn));
    b_mirror(1:pad1) = flip(b(1:pad1));
    B(:,t) = b_mirror;
end
end

