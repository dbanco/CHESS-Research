function [x, cgst] = solvemdbi_cg_OF(a, rho, b, tol, mit, isn,N,K,J,T,lambda2,U,V,Jterm)

% solvemdbi_cg_OF -- Solve a multiple diagonal block linear system with a
%                  scaled identity term using CG
%
%         The solution is obtained by independently solving a set of linear
%         systems of the form (see wohlberg-2016-efficient)
%s
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

if nargin < 6
  isn = [];
end
if nargin < 5 || isempty(mit)
  mit = 1000;
end
if nargin < 4 || isempty(tol)
  tol = 1e-5;
end

xsz = [1,N,K*J,T];
ah = conj(a);
Aop = @(u) sum(bsxfun(@times, a, u), 3);
Ahop = @(u) bsxfun(@times, ah, u);
AhAvop = @(u) vec(Ahop(Aop(reshape(u, xsz))));
AhAvop2 = @(u) vec(fft2(opticalFlowOp(ifft2(reshape(u,xsz),'symmetric'),U,V,K,1,Jterm) ));

wrn = warning('query','MATLAB:ignoreImagPart');
warning('off', 'MATLAB:ignoreImagPart');
if lambda2 == 0
    [xv,flg,rlr,pit,resvec] = pcg(@(u) AhAvop(u) + rho*u,...
    b(:), tol, mit, [], [], isn);
else
    [xv,flg,rlr,pit,resvec] = pcg(@(u) AhAvop(u) + rho*u + lambda2*AhAvop2(u),...
    b(:), tol, mit, [], [], isn);
end
warning(wrn.state, 'MATLAB:ignoreImagPart'); % 
cgst = struct('flg', flg, 'rlr', rlr, 'pit', pit);

x = reshape(xv, xsz);

return
