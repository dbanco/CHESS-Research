function [x, cgst] = solvemdbi_cg_smooth(ah, rho, b, tol, mit, isn, kappa)

% solvemdbi_ism -- Solve a multiple diagonal block linear system with a
%                  scaled identity term using CG
%
%         The solution is obtained by independently solving a set of linear
%         systems of the form (see wohlberg-2016-efficient)
%
%                  (rho I + a_0 a_0^H + a_1 a_1^H + ...) x = b
%
%         In this equation inner products and matrix products are taken along
%         the 3rd dimension of the corresponding multi-dimensional arrays; the
%         solutions are independent over the 1st and 2nd (and 4th, if
%         non-singleton) dimensions.
%
% Usage:
%       x = solvemdbi_cg(ah, rho, b);
%
% Input:
%       ah          Multi-dimensional array containing a^H
%       rho         Scalar rho
%       b           Multi-dimensional array containing b
%       tol         CG tolerance
%       mit         CG maximum iterations
%       isn         CG initial solution
%
% Output:
%       x           Multi-dimensional array containing linear system solution
%       cgst        CG status structure
%
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-05-14
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'License' file distributed with
% the library.

if nargin < 6,
  isn = [];
end
if nargin < 5 || isempty(mit),
  mit = 1000;
end
if nargin < 4 || isempty(tol),
  tol = 1e-5;
end


asz = [size(ah,1) size(ah,2) size(ah,3) 50];
a = conj(ah);
Aop = @(u) sum(bsxfun(@times, ah, u), 3);
Ahop = @(u) bsxfun(@times, a, u);
AhAvop = @(u) vec(Ahop(Aop(reshape(u, asz))));
DtDxop = @(u) DtDxfunc(u,asz);

wrn = warning('query','MATLAB:ignoreImagPart');
warning('off', 'MATLAB:ignoreImagPart');
[xv,flg,rlr,pit] = pcg(@(u) AhAvop(u)+rho*u+kappa*DtDxop(u), b(:), tol, mit, [], [], isn);
warning(wrn.state, 'MATLAB:ignoreImagPart');
cgst = struct('flg', flg, 'rlr', rlr, 'pit', pit);

x = reshape(xv, asz);

return

function y = Dxop(x)
[N1,N2,K,T] = size(x);
xt = reshape(x,[T,N1*N2*K]);
y = xt(1:(T-1),:) - xt(2:T,:);
return

function out = DtDxfunc(u,asz)
N1=asz(1);N2=asz(2);K=asz(3);T=asz(4);
x = reshape(u,asz);
xt = reshape(ifft2(x,'symmetric'),[T,N1*N2*K]);
y = xt(1:(T-1),:) - xt(2:T,:);
z = [y(1,:);...
    -y(1:(T-2),:)+y(2:(T-1),:);...
    -y(T-1,:)];
out = vec(fft2(reshape(z,asz),N1,N2));
return
