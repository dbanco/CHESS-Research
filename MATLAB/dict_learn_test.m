 %% Load data 
% 311 Ring dataset
dataset = fullfile('D:','CHESS_data','al7075_311_polar_reduced');

X = zeros(41*2048,185*5);
k = 1;
% load polar images
for load_step = 0:4
	for img = 0:184
        str1 = sprintf('%i',load_step);
        str2 = sprintf('%i',img);
        fileName = ['polar_image_',str1,'_',str2,'.mat'];
        fileDir = fullfile(dataset,fileName);
        load(fileDir)
        X(:,k) = polar_image(:);
        k = k + 1;
    end
end

%% Dictionary Learning

% initialize
Qh = [40,40,40,40,40]; Ph=[5,5,5,5,5];
Qv = [10,10,10,10,10]; Pv=[2,2,2,2,2];
sid = initSID2D(Qv,Qh,Pv,Ph,[]); 
plotSID2D(sid, 1);
D = makeSID2Dmatrix(sid,41,512); 
s = 200;
noIt = 500;
K = 10;

%% solve
for it = 1:noIt
    fprintf('Iter %i\n',it)
    W = sparseapprox(X(1:512*41,50:100:500), D, 'mexOMP', 'tnz',s);
    R = X(1:512*41,50:100:500) - D*W;
    for k=1:K
        fprintf('k: %i\n',k)
        I = find(W(k,:));
        Ri = R(:,I) + D(:,k)*W(k,I);
        [U,S,V] = svds(Ri,1,'L');
        % U is normalized
        D(:,k) = U;
        W(k,I) = S*V';
        R(:,I) = Ri - D(:,k)*W(k,I);
    end    
end