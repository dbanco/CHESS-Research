%% Load fit results
resultDir = 'D:\CHESS_data\al7075_311_alpha_beta_part2\*.mat';
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
xlabel('log(alpha)')
ylabel('log(# of coefs)')

figure(2)
plot(coefs,fitErr,'o') 
ylabel('Error')
xlabel('log(# of coefs')

figure(3)
semilogx(alphap(betap==min(betap)),fitErr(betap==min(betap)),'o')  
xlabel('log(alpha)')
ylabel('Error')

[task(betap==min(betap)) , betap(betap==min(betap))]
figure(4)
loglog(betap,coefs,'o')  
xlabel('log(beta)')
ylabel('log(# of coefs)')

figure(5)
semilogx(betap(alphap == 1),coefs(alphap == 1),'o')
xlabel('log(beta)')
ylabel('# of coefs')
