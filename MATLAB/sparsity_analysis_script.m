%% Script to plot fit sparsity and error
resultDir = 'D:\CHESS_data\al7075_311_alphap_weight_fit\*.mat';
files = dir(resultDir);

coefs = zeros(numel(files),1);
alphap = zeros(numel(files),1);
fitErr = zeros(numel(files),1);
for i = 1:numel(files)
    fprintf('%i of %i files \n',i,numel(files))
    load(fullfile(files(i).folder,files(i).name))
    coefs(i) = sum(x_hat(:)>0);
    alphap(i) = P.alphap;
    fitErr(i) = err(end);
end

figure(1)
loglog(alphap,coefs,'o')  
xlabel('log(alpha)')
ylabel('log(# of coefs)')

figure(2)
plot(coefs,fitErr,'o') 
ylabel('Error')
xlabel('log(# of coefs')

figure(3)
semilogx(alphap,fitErr,'o')  
xlabel('log(alpha)')
ylabel('Error')