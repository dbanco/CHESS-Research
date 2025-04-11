function [x, cgst] = solvemdbi_cgls_OF_gpu_zpad(a, rho, b1, b2, tol, mit, isn,N,M,K,J,T,lambda2,U,V)

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

xsz = [1,N+M-1,K*J,T];
% xszPad = [1,N+2*round(M/2),K*J,T];
a = gpuArray(complex(a));
ah = gpuArray(complex(conj(a)));


wrn = warning('query','MATLAB:ignoreImagPart');
warning('off', 'MATLAB:ignoreImagPart');
n=numel(isn);
if lambda2 == 0
    b = gpuArray(complex([b1(:);b2(:)]));
    m = numel(b);
    [xv,flg,rlr,pit] = cgls(@(ind,u,m,n) Aops1(ind,u,m,n,a,ah,xsz,M,rho),...
    0,b(:),m,n,mit,tol,1);
else
    b = gpuArray(complex([b1(:);b2(:);zeros((N+M-1)*K*J*T,1)]));
    m = numel(b);
    [xv,flg,rlr,pit] = cgls(@(ind,u,m,n) Aops2(ind,u,m,n,a,ah,xsz,M,rho,U,V,K,lambda2),...
    0,b(:),m,n,mit,tol,1);
end
warning(wrn.state, 'MATLAB:ignoreImagPart'); % 
cgst = struct('flg', flg, 'rlr', rlr, 'pit', pit);

x = reshape(xv, xsz);

end

function out = Aops1(ind,u,m,n,a,ah,xsz,M,rho)
    if ind == 1 
        u = reshape(u, xsz);
        out1 = sum(pagefun(@times, a, u), 3);
        out = [vec(out1);rho*vec(u)];
    elseif ind == 2 
        numelb = xsz(2)*xsz(4);
        numelx = xsz(2)*xsz(3)*xsz(4);
        bsz = [1,xsz(2),1,xsz(4)];
        u1 = reshape(u(1:numelb),bsz);
        u2 = reshape(u(numelb+1:numelb+numelx), xsz);
        out1 = pagefun(@times, ah, fft2(maskPad(ifft2(u1,'symmetric'),M)));
        out = vec(out1) + rho*vec(u2);
    elseif ind == 0
        out = 'Aops1';
    end
end
function out = Aops2(ind,u,m,n,a,ah,xsz,M,rho,U,V,K,lambda2)
    if ind == 1 
        u = reshape(u, xsz);
        out1 = sum(pagefun(@times, a, u), 3);
        out2 = vec(fft2(opticalFlowOp(real(ifft2(reshape(u,xsz),'symmetric')),U,V,K,0) ));
        out = [vec(out1);rho*vec(u);lambda2*vec(out2)];
    elseif ind == 2
        numelb = xsz(2)*xsz(4);
        numelx = xsz(2)*xsz(3)*xsz(4);
        bsz = [1,xsz(2),1,xsz(4)];
        u1 = reshape(u(1:numelb),bsz);
        u2 = reshape(u(numelb+1:numelb+numelx), xsz);
        u3 = reshape(u(numelb+numelx+1:end), xsz);
        out1 = pagefun(@times, ah, fft2(maskPad(ifft2(u1,'symmetric'),M)));
        out2 = vec(fft2(opticalFlowOp(real(ifft2(reshape(u3,xsz),'symmetric')),U,V,K,2) ));
        out = vec(out1) + rho*vec(u2) + lambda2*vec(out2);
    elseif ind == 0
        out = 'Aops2';
    end
end

function u = maskPad(u,M)
u(:,1:M-1,:,:) = 0;
end
