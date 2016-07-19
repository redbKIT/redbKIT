function [IDEIM_sparse,PHI_sparse] = MDEIM(U_sparse, m)
%MDEIM Matrix Discrete Empirical Interpolation Method
%
%   [IDEIM_SPARSE] = MDEIM(U_SPARSE) given a sparse basis matrix
%   U_SPARSE of size Nh^2 x M, returns a M x 1 vector IDEIM_SPARSE of indices 
%   from 1 to Nh^2. 

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

fprintf('\nStart MDEIM ... ');
timeMDEIM = tic;

M  = size(U_sparse,2);
Nh = sqrt(size(U_sparse,1));

if nargin < 2 || isempty(m)
      m = M;
end

if m > M
      m = M;
end

[Z, ~] = find(U_sparse);
Z      = unique(Z);
U      = full(U_sparse(Z,:));
clear U_sparse;

%% Start Index selection
[~, indx_max]     =  max(abs(U(:,1)));
IDEIM(1,1)        =  indx_max;
PHI               =  U(:,1);

for i = 2 : size(U,2)
    
      p             =  min(i-1,m);
      A             =  PHI(IDEIM, 1:p);  
      b             =  U(IDEIM, i);
      c             =  A\b;
      r             =  U(:,i) - PHI*c;
      
      % find worst approximated entry
      [~, indx_max] = max(abs(r));
      IDEIM(i,1)    = indx_max(1);
      
      if i <= m
            PHI(:,i)    =  U(:,i);
      end
end

IDEIM_sparse = Z(IDEIM);

[~, jP, vP] = find(PHI);
PHI_sparse  = sparse(repmat(Z,size(PHI,2),1),jP,vP,Nh^2,size(PHI,2));

timeMDEIMf = toc(timeMDEIM);
fprintf('done in %2.2e s\n',timeMDEIMf);

end
