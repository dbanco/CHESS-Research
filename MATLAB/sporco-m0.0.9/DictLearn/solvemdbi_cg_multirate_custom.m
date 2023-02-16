function [x, cgst] = solvemdbi_cg_multirate_custom(ah, rho, b, tol, mit, isn,scales,NormVals)

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
U = size(scales,2);
M = size(b,2);
K = size(b,3);
N2 = size(ah,2);
asz = [1 M K];
% Aop = @(in) ifft2(sum(bsxfun(@times,ah,fft2( reSampleNu(N2,in,c1,c2,Ufactors) )),3),'symmetric');
% Ahop = @(in) reSampleNuTrans2(M,ifft2(sum(bsxfun(@times, conj(ah), fft2(in)), 4),'symmetric'),c1,c2,Ufactors);
AhAvop = @(in) vec(wrapAhA(reshape(in, asz),N2,scales,ah,M,NormVals));

% Aop = @(u) ifft2(sum(bsxfun(@times, ah, fft2(u)), 3),'symmetric');
% Ahop = @(u) ifft2(sum(bsxfun(@times, conj(ah), fft2(u)), 4),'symmetric');
% AhAvop = @(u) vec(Ahop(Aop(reshape(u, asz))));

wrn = warning('query','MATLAB:ignoreImagPart');
warning('off', 'MATLAB:ignoreImagPart');
[xv,flg,rlr,pit] = pcg(@(u) AhAvop(u)+rho*u, b(:), tol, mit, [], [], isn);
warning(wrn.state, 'MATLAB:ignoreImagPart');
cgst = struct('flg', flg, 'rlr', rlr, 'pit', pit);

x = reshape(xv, asz);

end

function out = wrapAhA(in,N2,scales,ah,M,NormVals)
A1 = reSampleCustomArray(N2,in,scales,NormVals);
Ain = sum(bsxfun(@times,ah,fft2( A1 )),3);
out = reSampleTransCustomArray(M,ifft2(sum(bsxfun(@times, conj(ah), Ain), 4),'symmetric'),scales,NormVals);
end
