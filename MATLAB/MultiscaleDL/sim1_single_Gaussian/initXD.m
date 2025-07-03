function opt = initXD(opt,N,M,K,J,T,Xtrue,Dtrue)
KJ = K*J;

if nargin <= 6
    Xtrue = 0;
    Dtrue = 0;
end

% Init Coefficents
switch opt.coefInit
    case 'zeros'
        opt.Y0 = zeros(1,N+M-1,KJ,T);
    case 'rand'
        opt.Y0 = rand(1,N+M-1,KJ,T);
    case 'true'
        opt.Y0 = Xtrue;

end

% Init Dictionary
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
switch opt.dictInit
    case 'flat'
        D0 = zeros(1,M,K);
        for k = 1:K
            i1 = round(M*(0.5 - 0.5*k/K) + 1);
            i2 = round(M*(0.5 + 0.5*k/K));
            D0(1,i1:i2,k) = k;
            D0(1,:,k) = Pnrm(D0(1,:,k));
        end
    case 'zeros'
        D0 = zeros(1,M,K);
    case 'rand'
        D0 = rand(1,M,K);
    case 'true'
        D0 = Dtrue;
end
opt.G0 = D0;