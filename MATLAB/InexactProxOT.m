function [Wd, r, c, Gamma, a_t] = InexactProxOT(r, c, D, beta, L)
%sinkhornKnoppTransport

% Remove 0 probability points from each distribution, normalize
Indr = r > 0;
r = r(Indr);
r = r./sum(r(:));

Indc = c > 0;
c = c(Indc);
c = c./sum(c(:));

D = D(Indr,Indc);
G = exp(-D/beta);
Gamma = ones(size(D));

b = ones(size(r))/numel(r);
a_t = {};
for i = 1:1000
    Q = G.*Gamma;
    for t = 1:L
       a = r./(Q*b);
       b = c./(Q'*a);
    end
    a_t{i} = a;
    Gamma_new = diag(a)*Q*diag(b);
    check = norm(Gamma_new-Gamma);
    if sum(abs(check ))<1e-8
        break
    end
    Gamma = Gamma_new;
end

Wd = sum(Gamma(:).*D(:));

end
