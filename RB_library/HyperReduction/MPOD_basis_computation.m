function [U_sparse, V, S] = MPOD_basis_computation(X_sparse, tolPOD, flag_svd)
%MPOD_basis_computation implements POD for vectorized sparse matrices
%
%   [U_SPARSE, S, V] = MPOD(X_SPARSE, TOLPOD) given a sparse snasphots 
%   matrix X_SPARSE of size Nh^2 x Ns and a tolerance (or number of POD modes 
%   to extract), returns a sparse POD basis U_sparse of size Nh^2 x M, 
%   the corresponding singular values S and a unitary matrix V of size M x
%   M.
%
%   [U_SPARSE, S, V] = MPOD(X_SPARSE, TOLPOD, FLAG_SVD) if FLAG_SVD = 1, 
%   employs a more accurate (but slower) thin SVD for the computation of 
%   the POD basis.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 3 || isempty(flag_svd)
    flag_svd = 1;
end

fprintf('\nStart MPOD ... ');
timeMPOD = tic;

Nh = sqrt(size(X_sparse,1));
Ns = size(X_sparse,2);

% check if the number of nonzeros changes over the columns 
[Z, ~] = find(X_sparse);
Z      = unique(Z);
    
X = full(X_sparse(Z,:));
clear X_sparse;

%% Compute eigen- or singular values decomposition of X
if flag_svd
    % Thin SVD
    [U, S, V] = svd(X,'econ');
else
    % method of snapshots: find eigenvalues of the correlation matrix
    C           = X' * X;
    [V,S,~]     = svd(full(C));    
    S           = sqrt(S);
end

S  = diag(S);

%% Extract POD modes up to tolPOD tolerancce
if tolPOD < 1
      S_N  = cumsum(S.^2);
      M    = find(tolPOD.^2>=1-S_N/sum(S.^2),1,'first');
else
      M    = min(tolPOD,size(V,2)); 
end
fprintf('%d modes selected,',M);

if flag_svd
    U   =  U(:, 1 : M);
    V   =  V(:, 1 : M);
else
    % Form dense POD basis
    V   =  V(:, 1 : M);
    U   =  X*V;
    for cc = 1 : M
        U(:,cc) = 1/(S(cc)) * U(:,cc);
    end
end

%% Form sparse POD basis
[~, jU, vU] = find(U);
U_sparse    = sparse(repmat(Z,M,1), jU, vU, Nh^2, M);

timeMPOD2   = toc(timeMPOD);
fprintf(' done in %2.2e s\n',timeMPOD2);

end
