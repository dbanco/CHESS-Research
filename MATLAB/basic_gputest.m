%% Basic GPU test script
gpuDevice

A1 = rand(3000,3000);
tic;
B1 = fft(A1);
time1 = toc;

A2 = gpuArray(A1);
tic;
B2 = fft(A2);
time2 = toc;