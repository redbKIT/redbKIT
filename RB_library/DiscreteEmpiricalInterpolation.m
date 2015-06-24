function [IDEIM, PHI, P] = DiscreteEmpiricalInterpolation(U, m)
%DISCRETEEMPIRICALINTERPOLATION performs DEIM algorithm
%
%   [IDEIM, PHI, P] = DISCRETEEMPIRICALINTERPOLATION(U) given a matrix U 
%   of M column snapshot vectorsm, returns: (1) a vector IDEIM of length M 
%   containing the indices selected by DEIM, (2) a basis PHI = U, (3) a 
%   restriction (mask) matrix P.
%
%   [IDEIM, PHI, P] = DISCRETEEMPIRICALINTERPOLATION(U, m) if m<M, performs 
%   a least-squares (gappy POD) interpolation by using only m out of M basis 
%   vectors. In this case, PHI = U(:,1:m).
%
%   Reference: NONLINEAR MODEL REDUCTION VIA DISCRETE EMPIRICAL
%   INTERPOLATION, S. Chaturantabut & D. Sorensen, SIAM J. SCI. COMPUT.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%
%   The original DEIM algorithm has been provided by 
%   David Amsallem (Standford University), January 2014.


M = size(U,2);

if nargin < 2 || isempty(m)
      m = M;
end

if m > M
      m = M;
end

n                 = size(U,1);

P                 = [];

[~, indx_max]     =  max(abs(U(:,1)));

IDEIM(1,1)        =  indx_max;
P                 =  [P,unit(IDEIM(1,1),n)];

PHI               =  U(:,1);

for i = 2 : size(U,2)
    
      p             =  min(i-1,m);
      
      A             =  access_matrix(PHI, IDEIM, 1:p);  
 
      B             =  access_matrix(U, IDEIM, i);
      
      c             =  A\B;

      r             =  U(:,i) - linear_combo(PHI,c);
            
      [~, indx_max] = max(abs(r));
      
      IDEIM(i,1)    = indx_max(1);
      P             = [P,unit(IDEIM(i,1),n)];
      
      if i <= m
            PHI(:,i)    =  U(:,i);
      end
end



end

function e = unit(i,n)

e = sparse(n,1);
e(i,1) = 1;

end

function [C] = dot_matrix(A,B)

[~, M1] = size(A);
[~, M2] = size(B);

C      = sparse(M1,M2);
for i = 1 : M1
      for j = 1 : M2
            C(i,j) = dot(A(:,i), B(:,j));
      end
end

end

function [u] = linear_combo(PHI,c)

[n,  M]  = size(PHI);
[M2,tmp] = size(c);

u      = sparse(n,1);
for j = 1 : M
    u = u + PHI(:,j)*c(j);
end

end

function C = access_matrix(A, I, J)

tmp = A(:,J);
C   = tmp(I,:);

end

