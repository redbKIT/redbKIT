function [V,Sigma,PSI] = POD_basis_computation(u, Xnorm, N_tol, D)
%POD_BASIS_COMPUTATION
%
% INPUT: 
%         u: snapshots matrix
%         Xnorm: matrix norm
%         N_tol: either number of basis to extract or some tolerance on
%                 the energy to capture
%         D: quadrature weights norm
%
% OUTPUT: 
%         V      : reduced basis matrix
%         PSI    : selected eigenvectors
%         Sigma  : singular values

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

n  = size(u,1);
ns = size(u,2);

if nargin < 4 || isempty(D)
    D = speye(ns,ns);
end

H_D   = sqrt(D);
u_D   = full( u*H_D );

% Form Correlation matrix
if nargin < 2 || isempty(Xnorm)
      C = u_D'*u_D;        % u'*u;
else
      C = u_D'*(Xnorm*u_D);% u'*Xnorm*u;
end

[PSI, Sigma]  = svd(full(C));
Sigma         = diag(Sigma);     % squared singular values

% Extract Modes
if N_tol < 1
      
      Sigma_N  = cumsum(Sigma);
  
      N        = find(N_tol.^2>=1-Sigma_N/sum(Sigma),1,'first');
else
      
      N        = min(N_tol,size(PSI,2));
    
end

fprintf('\n   POD_basis_computation: %d modes selected \n',N);

PSI   =  PSI(:,1:N);

% V = [];
% for cc = 1 : size(PSI,2)    
%     V = [V 1/sqrt(Sigma(cc))*linear_combo(u_D,PSI(:,cc))];
% end

% Form POD basis
V   =  u_D*PSI;
for cc = 1 : N
    V(:,cc) = 1/(Sigma(cc)) * V(:,cc);
end

Sigma = sqrt(Sigma);

end


function [u] = linear_combo(PHI,c)

[n,  M]  = size(PHI);
[M2,tmp] = size(c);

u      = sparse(n,1);
for j = 1 : M
      u = u + PHI(:,j)*c(j);
end

end
