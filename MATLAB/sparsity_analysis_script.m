%% Load fit results
resultDir = 'D:\CHESS_data\al7075_polar_fit_ab_9_23\*.mat';
files = dir(resultDir);

coefs = zeros(numel(files),1);
alphap = zeros(numel(files),1);
betap = zeros(numel(files),1);
fitErr = zeros(numel(files),1);
task = zeros(numel(files),1);
for i = 1:numel(files)
    fprintf('%i of %i files \n',i,numel(files))
    load(fullfile(files(i).folder,files(i).name))
    coefs(i) = sum(x_hat(:)>0);
    alphap(i) = P.alphap;
    betap(i) = P.betap;
    task(i) = P.task;
    fitErr(i) = err(end);
end
%% Plot figures
figure(1)
loglog(alphap,coefs,'o')  
xlabel('alpha')
ylabel('# of coefs')

figure(2)
plot(coefs,fitErr,'o') 
ylabel('Error')
xlabel('# of coefs')

figure(3)
semilogx(alphap(betap==min(betap)),fitErr(betap==min(betap)),'o')  
xlabel('alpha')
ylabel('Error')

figure(33)
semilogx(alphap,fitErr,'o')  
xlabel('alpha')
ylabel('Error')

figure(44)
semilogx(betap,fitErr,'o')  
xlabel('beap')
ylabel('Error')

figure(4)
loglog(betap,coefs,'o')  
xlabel('beta')
ylabel('# of coefs')

figure(5)
idx = 78;
plot(betap(alphap == alphap(idx)),coefs(alphap == alphap(idx)),'o')
xlabel('beta')
ylabel('# of coefs')
