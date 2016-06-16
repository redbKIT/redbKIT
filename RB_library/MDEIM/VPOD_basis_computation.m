function [U, V, S] = VPOD_basis_computation(X, Xnorm, tolPOD, flag_svd)
%VPOD_BASIS_COMPUTATION implements POD for vectors

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 4 || isempty(flag_svd)
    flag_svd = 1;
end

fprintf('\nStart VPOD ... ');
timeVPOD = tic;

X        = full(X);

%% Compute eigen- or singular values decomposition of X
if flag_svd
    % Thin SVD
    if nargin < 2 || isempty(Xnorm)
        [U, S, V]   = svd(X,'econ');
    else
        [L, flag, perm] = chol(Xnorm,'lower','vector');
        X               = L'*X(perm,:);
        invp            = 0*perm;
        invp(perm)      = 1:length(perm);
        X               = X(invp,:);
        [U, S, V]       = svd(X,'econ');
        U               = L' \ U(perm,:);
        U               = U(invp,:);
    end
    
else
    % method of snapshots: find eigenvalues of the correlation matrix
    if nargin < 2 || isempty(Xnorm)
        C =  ( X'* X ) ;
    else
        C =  ( X'* (Xnorm * X) );
    end
    [V, S, ~]     = svd(full(C));    
    S             = sqrt(S);
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
    % Form POD basis
    V   =  V(:, 1 : M);
    U   =  X*V;
    for cc = 1 : M
        U(:,cc) = 1/(S(cc)) * U(:,cc);
    end
end

timeVPOD2 = toc(timeVPOD);
fprintf(' done in %2.2e s\n',timeVPOD2);

end


