% Analyze results, should I run again?
c1s = [results.c1];
c2s = [results.c2];
outputs = [results.output];
dictErr = [results.dictErr];


inFracs = c1s./c2s;
outFracs = outputs(1,:)./outputs(2,:);

figure(1)
plot(inFracs,outFracs,'-o')
hold on
plot(inFracs,inFracs,'-')

figure(2)
plot(abs(inFracs-outFracs),'-o')
