
% Composite op: D3 * U2, Adjoint: U3 * D2
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);
c = scales{1};



k = 1;
count = 0;
M = 71;
for N = 117:2:121   
    % Create random signals
    d = zeros(1,M,1);
    center = ceil((M+1)/2);
    d(1,16:56,1) = [1:20,21,20:-1:1];
    
    x = randn(1,N+M-1,K*J);
    xh = conj(x);

    
    
    [Md,NormVals,Shifts] = reSampleCustomArrayCenter3(N,d,scales,center); 
    y = rand(1,N,K*J);
    Mdpad = padarray(Md,[0 M-1 0 0],0,'post');
    convMd = ifft2(sum(bsxfun(@times,x,fft2( Mdpad )),3),'symmetric');
    convMd(:,1:M-1,:,:) = 0;
    % convtconvMd = ifft2(sum(bsxfun(@times, xh, fft2(convMd)), 4),'symmetric');
    % MtconvtconvMd = reSampleTransCustomArrayCenter(M,convtconvMd,scales,center,NormVals);
    Ad = convMd;
   
    ypad = padarray(y,[0 M-1 0 0],0,'post');
    convty = ifft2(sum(bsxfun(@times, xh, fft2(ypad)), 4),'symmetric');
    Mtconvty = reSampleTransCustomArrayCenter3(M,convty,scales,center,NormVals,Shifts);
    % convconvtMty = ifft2(sum(bsxfun(@times,x,fft2( convtMty )),3),'symmetric');
    Aty = Mtconvty;

    % === TEST ===
    lhs = sum(Ad .* ypad,'all');
    rhs = sum(d .* Aty,'all');

    % fprintf('Composite adjoint test:\n');
    % fprintf('lhs = %.12e, rhs = %.12e, error = %.2e\n', ...
    %     lhs, rhs, abs(lhs - rhs));
    
    if abs(lhs - rhs) > 1e-10
        fprintf('N = %i \n',N)
        count = count + 1;
    end
    k = k + 1;
end
fprintf('%i/%i errors \n',count,k-1)
