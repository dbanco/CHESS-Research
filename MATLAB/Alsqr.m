function out = Alsqr(in,type,U,xf,N,T,K,sigma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
if type == 1
    in1 = reshape(in,[1,N,K,1]);
    out2 = ifft2(sum(bsxfun(@times,xf,fft2(Pnrm(decimate(in1,U)))),3),'symmetric');
    out = [out2(:); sqrt(sigma)*in];
elseif type == 2
    y = reshape(in(1:N*T),[1,N,1,T]);
    AXy = decimateTrans(ifft2(sum(bsxfun(@times, conj(xf), fft2(y)), 4),'symmetric'),U);
    out = AXy(:) + sqrt(sigma)*in(N*T+1:end);
else
    error('Valid values for type are 1 and 2')
end

