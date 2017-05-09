%% Basic GPU test script
tic;
A4 = rand(3000,3000);
B4 = fft(A4);
time4 = toc;

tic;
A5 = rand(3000,3000,'gpuArray');
B5 = fft(A5);
B5 = gather(B5);
time5 = toc;

speedUp = time4/time5;
disp(speedUp);