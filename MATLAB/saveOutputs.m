function saveOutputs(topDir,ending,Y,U,D1,D2,X1,X2,Vals,param)
%% Dir to save figures
save(fullfile(topDir,['variables',ending,'.mat']),'U','D1','D2','X1','X2','Vals','param')

%% Visualize factors in time order
% Sort by time occurence
[T,R]= size(U{3});
timeW = 1:T;
timeCoM = zeros(R,1);
for i = 1:R
    timeCoM(i) = timeW*U{3}(:,i)/sum(U{3}(:,i),'all');
end
[~,ind] = sort(timeCoM);

f1 = figure; % View angle components
for i = 1:R
    subplot(5,4,i)
    j = ind(i);
    angComp = U{1}(:,j)*U{2}(:,j)';
    imagesc(angComp)
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    title(['t= ',num2str(timeCoM(j))])
end
saveas(f1,fullfile(topDir,['UU12',ending,'.png']))

U3sort = U{3};
for i = 1:R
    j = ind(i);
    U3sort(:,i) = U{3}(:,j);
end
f2 = figure;
imagesc((U3sort'))
saveas(f2,fullfile(topDir,['U3',ending,'.png']))

% f3 = figure;
% for i = 1:R
%     subplot(5,4,i)
%     j = ind(i);
%     k = floor(timeCoM(j));
%     imagesc(Y(:,:,k))
%     title(['t= ',num2str(k)])
% end
% saveas(f3,fullfile(topDir,['DataEx',ending,'.png']))

% View learned dictionaries
f4 = plotDictionary(D1);
saveas(f4,fullfile(topDir,['D1',ending,'.png']))
f5 = plotDictionary(D2);
saveas(f5,fullfile(topDir,['D2',ending,'.png']))

Y1 = recon(D1,X1);
Y2 = recon(D2,X2);

% Error in R
rel_err = zeros(R,2);
sparsity = zeros(R,2);
dataNorms = zeros(R,1);
% vdfs1 = squeeze(squeeze(sum(sum(X1,1),2)));
% vdfs2 = squeeze(squeeze(sum(sum(X2,1),2)));
for j = 1:R
    i = ind(j);
   rel_err(j,1) = norm(Y1(:,i)-U{1}(:,i))/norm(U{1}(:,i));
   rel_err(j,2) = norm(Y2(:,i)-U{2}(:,i))/norm(U{2}(:,i));
   sparsity(j,1) = sum(X1(:,:,:,i)>0,'all');
   sparsity(j,2) = sum(X2(:,:,:,i)>0,'all');
%    vdfs1(:,j) = vdfs1(:,i)/max(vdfs1(:,i));
%    vdfs2(:,j) = vdfs2(:,i)/max(vdfs2(:,i));
end

f6 = figure;
subplot(2,1,1)
plot(rel_err)
legend('D1X1-U1','D2X2-U2')
title('Relative Error norm(DiXi-Ui')

subplot(2,1,2)
plot(sparsity)
legend('X1','X2')
title('Number of Nonzeros in Xi')

% subplot(4,1,3)
% imagesc(vdfs1)
% title('2\theta VDFs')
% xlabel('R')
% 
% subplot(4,1,4)
% imagesc(vdfs2)
% title('\eta VDFs')
% xlabel('R')

saveas(f6,fullfile(topDir,['err_l1',ending,'.png']))

% View CP factors
f7 = plotDictionary(U{1}(:,ind));
saveas(f7,fullfile(topDir,['U1',ending,'.png']))
f8 = plotDictionary(U{2}(:,ind));
saveas(f8,fullfile(topDir,['U2',ending,'.png']))

% View objective
f9 = figure;
valNames = {'Itn','Obj','ncpErr','D1Err',...
          'D2Err','L3Err','D1L1','D2L1'};
for i = 2:8
subplot(4,2,i-1)
keep = sum(Vals,2)>0;
semilogy(Vals(keep,i))
title(valNames(i))
end
saveas(f9,fullfile(topDir,['obj',ending,'.png']))


% View reconstructions
reconDir = fullfile(topDir,['recon',ending]);
mkdir(reconDir)
Yhat = double(full(ktensor(U)));
N2 = size(Y,2);
f10 = figure;
for t = 1:T
    maxY = max(Y(:,:,t),[],'all');
    imagesc([Y(:,:,t); maxY*ones(1,N2); Yhat(:,:,t)])
    saveas(f10,fullfile(reconDir,['recon_',num2str(t),'.png']))
    relErr = norm(Y(:,:,t)-Yhat(:,:,t),'fro')/norm(Y(:,:,t),'fro');
    title(sprintf('Relative error: %2.2f',relErr))
end
close all
end

function Y = recon(D,X)
Df = fft2(D,1,size(X,2));
Xf = fft2(X);
Yf = sum(bsxfun(@times,Df,Xf),3);
Y = squeeze(real(ifft2(Yf)));
end