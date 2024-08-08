function y = lowpassM(x,m,trans)
N = numel(x);

% Convolution kernel (linear interpolation)
if m > 1
%     f = zeros(1,N);
%     f(1) = 1;
%     f(end-m+2:end) = linspace(1/m,1-1/m,m-1);
%     f(2:m) = linspace(1-1/m,1/m,m-1);
    f = [linspace(1/m,1-1/m,m-1),1,linspace(1-1/m,1/m,m-1)];
    f = f/norm(f);
    if isgpuarray(x)
        f = gpuArray(f);
    end
%     y = real(ifft( fft(complex(x)).*fft(f)));
    y = convn(x,f,'same');
else
    y = x;
end

end