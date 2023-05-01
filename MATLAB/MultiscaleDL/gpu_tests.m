A = randn(100,300,40,100);
dictA = randn(100,300,40);
B = gpuArray(complex(randn(1,300,40,100)));
dictB = gpuArray(randn(1,300,40));

tic
for t=1:10
    Af = fft2(A);
    Daf = fft2(dictA);
    outA = ifft2(bsxfun(@times,Af,Daf),'symmetric');
end
toc

tic
for t=1:10
    Bf = fft2(B);
    Dbf = fft2(dictB);
    outB = ifft2(pagefun(@times,Bf,Dbf),'symmetric');
end
toc