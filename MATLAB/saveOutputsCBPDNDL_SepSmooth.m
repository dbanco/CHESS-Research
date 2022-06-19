function saveOutputsCBPDNDL_SepSmooth(topDir,ending,Y,D,X,Dx,Dy,optinf,lambda)
%% Dir to save figures
save(fullfile(topDir,['variables',ending,'.mat']),'D','X','Dx','Dy','optinf','lambda')

% Sort by time usage
% Sort by time occurence
[~,~,R,T]= size(X);
vdf_time = squeeze(sum(sum(X,1),2));

timeW = 1:T;
timeCoM = zeros(R,1);
for r = 1:R
    timeCoM(r) = timeW*vdf_time(r,:)'/sum(vdf_time(r,:),'all');
end
[~,ind] = sort(timeCoM);

f1 = figure;
Dsort = zeros(size(D));
for r = 1:R
    Dsort(:,:,r) = D(:,:,ind(r));
end

imshow(tiledict(Dsort),'Border','tight'); 
colormap jet
f1.Position = [100 100 540 400];
saveas(f1,fullfile(topDir,['dictionary',ending,'.png']))

% Plot functional value evolution
f2 = figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');
saveas(f2,fullfile(topDir,['obj',ending,'.png']))

% View reconstructions
reconDir = fullfile(topDir,['recon',ending]);
mkdir(reconDir)
Yhat = recon(D,X);
[N1,N2,R,T] = size(X);
f3 = figure;
T = size(Y,3);
Xplace = zeros(N1,N2,T); 
for t = 1:T   
    maxY = max(Y(:,:,t),[],'all');
    maxYhat = max(Yhat(:,:,t),[],'all');
    imshow([Y(:,:,t)/maxY; ones(1,N2);...
            Yhat(:,:,t)/maxYhat;...
            ],'Border','tight')
    
    relErr = norm(Y(:,:,t)-Yhat(:,:,t),'fro')/norm(Y(:,:,t),'fro');
    title(sprintf('Relative error: %2.2f',relErr))
    colormap jet
    f3.Position = [100 100 540 400];
    saveas(f3,fullfile(reconDir,['recon_',num2str(t),'.png']))
end

% Dictionary entry distribution
f4 = figure;
vdf_time_sort = zeros(size(vdf_time));
for r = 1:R
    vdf_time_sort(r,:) = vdf_time(ind(r),:);
end
imshow(vdf_time_sort,'Border','tight', 'InitialMagnification', 500)
colormap jet
xlabel('time')
ylabel('Dict Entry #')
title('Dictionary Entry usage in Time')
f4.Position = [100 100 540 400];
saveas(f4,fullfile(topDir,['vdf_time',ending,'.png']))

% Unfolded in time to see smoothness
f5 = figure;
Xtime = reshape(X,[T,N1*N2*R]);
imagesc(log(Xtime))
xlabel('N1*N2*K')
ylabel('Time')
f5.Position = [100 100 1000 500];
saveas(f5,fullfile(topDir,['X_unfold',ending,'.png']))


% Placement of dictionary entries
% f5 = figure;
% N2 = size(X,2);
% for t = 25
%     showCoef = [];
%     coefs = X(:,:,:,t);
%     maxCoef = max(coefs(:));
%     for r = 1:8
%         showCoef = [showCoef;X(:,:,r,t)/maxCoef;ones(1,N2)];
%     end
%     imshow(showCoef/maxCoef,'Border','tight', 'InitialMagnification', 500)
%     colormap jet
% end

close all
end


function Y = recon(D,X)
Df = fft2(D,size(X,1),size(X,2));
Xf = fft2(X);
Yf = sum(bsxfun(@times,Df,Xf),3);
Y = squeeze(real(ifft2(Yf)));
end
