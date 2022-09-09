function PhifGf = scaleDict(Phif,Gf)
[N1,N2,U] = size(Phif);
K = size(Gf,3);
PhifGf = zeros(N1,N2,U*K);
for k = 1:K
    for u = 1:U
        PhifGf(:,:,(u+(k-1)*U)) = bsxfun(@times,Phif(:,:,u),Gf(:,:,k));
    end
end
PhiG = ifft2(PhifGf,'symmetric');
for i = 1:(U*K)
    PhiG(:,:,i) = PhiG(:,:,i)/norm(PhiG(:,:,i));
end
PhifGf = fft2(PhiG);
end

