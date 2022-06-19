function y = lwpass4(x,trans)
N = numel(x);
% f = exp(-1)*ones(1,4);
f = exp(-1)*[0.5 1 0.5]
M = numel(f);

f4 = zeros(1,N);
if trans
    f4(1) = f(1);
    f4(N-(M-2):N) = f(M:-1:2);
else
    f4(1:M) = f;
end
y = ifft( fft(x).*fft(f4));
end

