function [IEIM, PHI, muEIM] = EmpiricalInterpolation(fun, Omega_h, Xi_train, tol, Max_Iter)
%EMPIRICALINTERPOLATION performs EIM algorithm
%
%   [IEIM, PHI, MUEIM] = EMPIRICALINTERPOLATION(FUN, OMEGA_H, XI_TRAIN, TOL) perform EIM on
%   the function FUN. FUN can be an m-function, an inline function or a
%   function handle, e.g. fun = @(x,y)(x.^y). OMEGA_H is a matrix of dimension
%   NumberOfXPoints x XDimension, XI_TRAIN is a matrix of dimension
%   NumberOfYPoints x YDimension. TOL is the desired tolerance for the
%   stopping criterion. The function returns: (1) a vector IEIM of length M
%   containing the indices selected by EIM, (2) a basis PHI of dimension
%   NumberOfXPoints x M, (3) the selected parameters values MUEIM.
%
%   [IEIM, PHI, MUEIM] = EMPIRICALINTERPOLATION(FUN, OMEGA_H, XI_TRAIN, TOL, MAX_ITER)
%   at most MAX_ITER iterations are performed
%
%   Reference: AN 'EMPIRICAL INTERPOLATION METHOD: APPLICATION TO
%   EFFICIENT REDUCED-BASIS DISCRETIZATION OF PARTIAL DIFFERENTIAL
%   EQUATIONS, M. Barrault, Y. Maday, N.C. Nguyen, A.T. Patera

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

fprintf('\n** Start Empirical Interpolation **\n')

Ns   = size(Xi_train, 1);
Nq   = size(Omega_h, 1);

if nargin < 5 || isempty(Max_Iter)
    Max_Iter = Ns;
end

M           = 0;
EIM_error   = tol + 1;

% find first basis function
e_M = zeros(Ns,1);
for j = 1 : Ns
    e_M(j) = norm( fun(Omega_h, Xi_train(j,:)) , Inf);
end
[~ , new_indx_basis] = max(e_M);
muEIM(1,:)           = Xi_train(new_indx_basis,:);              

r   =  fun( Omega_h, muEIM(1,:) );

while EIM_error > tol && M < Max_Iter
    
    M = M + 1;
    
    % find M-th interpolation point x^M  = x(IEIM(M))
    [~, indx_max]  = max(abs(r));
    IEIM(M)        = indx_max(1);
    PHI(:,M)       = r / r(IEIM(M));
    
    % find (M+1)-th basis function
    B      =  PHI(IEIM, 1:M);
    e_M    = zeros(Ns,1);
    for j = 1 : Ns
        
        fun_j    = fun(Omega_h, Xi_train(j,:));
        
        r        = fun_j - PHI* (B \ fun_j(IEIM));
        
        e_M(j)   = norm(r, Inf)/norm(fun_j, Inf);
        
    end
    
    [EIM_error, new_indx_basis] = max(e_M);
    
    fprintf('iter = %d, max e_M = %e\n', M, EIM_error);
    
    % compute residual
    muEIM(M+1,:) = Xi_train(new_indx_basis,:);      
    fun_next     = fun( Omega_h, muEIM(M+1,:) );
    r            = fun_next - PHI * (B \ fun_next(IEIM));
     
end

fprintf('\n** End Empirical Interpolation **\n')

end
