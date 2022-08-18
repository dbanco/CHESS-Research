function [x_hat] = IOCG( A0ft_stack,b,yk,vk,params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tau = 1e-4;
lambda = params.lambda;
beta = params.beta;
maxIter = params.maxIter;

[N,M] = size(A0ft_stack);
b_ft = fft(b);
xk = AtR_ft_1D(A0ft_stack,b);
dk = ones(N,1);
s = A0ft_stack./repmat(b_ft,[1,M]); 

for i = 1:maxIter
    yk = inner(xk,dk);
    [xk,dk] = project(yk);   
    stopCondition = (min(yk) < tau);
    if stopCondition
        break;
    end
end

x_hat = xk;

function [xkp1,dkp1] = project(yk)
dkp1 = yk >= 0;
xkp1 = yk.*dkp1;
end

function [yk,kin] = inner(x0,d)

xk = fft(x0);
b_ft = fft(b);
ax_ft = Ax_1D(A0ft_stack,xk);
rk = b_ft - ax_ft;
pk = AtR_1D(A0ft_stack,rk,d);
qk = pk;

t = s.*xk;
deltakm1 = N^2*sum(rk(:).^2)./(N^2 - sum(t(:)))^2;

for ii = 1:100
    zk = Ax_1D(A0ft_stack,pk);
    alphak = sum(qk(:).^2)/sum(zk(:).^2);
    xkp1 = xk + alphak*pk;
    rkp1 = rk - alphak*zk;
    qkp1 = AtR_1D(A0ft_stack,rkp1,d);
    betak = sum(qkp1(:).^2)/sum(qk(:).^2);
    pkp1 = qkp1 + betak*pk;
    [chi,deltakm1] = isc(xkp1,rkp1,deltakm1);
    if chi
        break
    end
    xk = xkp1;
    rk = rkp1;
    qk = qkp1;
    pk = pkp1;
end

yk = ifft(xk);
end

function [chi,deltak] = isc(xkp1,rkp1,deltakm1)
    t = s.*xkp1;
    deltak = sum(rkp1(:).^2)./(1 - sum(t(:)))^2;
    chi = (deltak < deltakm1);       
end

end

function AtR = AtR_1D(A,R,d)
    AtR = zeros(size(A));
    for j = 1:size(A,2)
        y = A(:,j).*R(:);
        AtR(:,j) = shift1D(y,2);
    end
    AtR = fft(d.*ifft(AtR));
end

function Ax = Ax_1D(A,x)
    Ax = zeros(size(A,1),1);
    for k = 1:size(x,2)
            Ax = Ax + A(:,k).*x(:,k);
    end
end