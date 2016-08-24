function [IDEIM, P, PHI] = DEIM(U, m)
%DISCRETEEMPIRICALINTERPOLATION performs DEIM algorithm
%
%   [IDEIM, P, PHI] = DISCRETEEMPIRICALINTERPOLATION(U) given a matrix U 
%   of M column snapshot vectorsm, returns: (1) a vector IDEIM of length M 
%   containing the indices selected by DEIM, (2) a basis PHI = U, (3) a 
%   restriction (mask) matrix P.
%
%   [IDEIM, P, PHI] = DISCRETEEMPIRICALINTERPOLATION(U, m) if m<M, performs 
%   a least-squares (gappy POD) interpolation by using only m out of M basis 
%   vectors. In this case, PHI = U(:,1:m).
%
%   Reference: NONLINEAR MODEL REDUCTION VIA DISCRETE EMPIRICAL
%   INTERPOLATION, S. Chaturantabut & D. Sorensen, SIAM J. SCI. COMPUT.

%   The original DEIM algorithm has been provided by 
%   David Amsallem (Standford University), January 2014.

fprintf('\nStart DEIM ... ');

timeDEIM = tic;

Nh = size(U,1);
M  = size(U,2);

if nargin < 2 || isempty(m)
      m = M;
end

if m > M
      m = M;
end

%% Start Index selection
[~, indx_max]     =  max(abs(U(:,1)));
IDEIM(1,1)        =  indx_max;

PHI               =  U(:,1);
P                 =  [];
P                 =  [P,unit(IDEIM(1,1),Nh)];

for i = 2 : size(U,2)
    
      p             =  min(i-1,m);
      
      A             =  PHI(IDEIM, 1:p);  
 
      B             =  U(IDEIM, i);
      
      c             =  A\B;

      r             =  U(:,i) - PHI*c;
            
      [~, indx_max] = max(abs(r));
      
      IDEIM(i,1)    = indx_max(1);
      P             = [P,unit(IDEIM(i,1),Nh)];
      
      if i <= m
            PHI(:,i)    =  U(:,i);
      end
end

timeDEIMf = toc(timeDEIM);

fprintf('done in %2.2e s\n',timeDEIMf);

end
 
function e = unit(i,n)

e = sparse(n,1);
e(i,1) = 1;

end
