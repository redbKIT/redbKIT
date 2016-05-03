function [uN, uNh, betaN] = solve_RBsystem(ROM, mu, N, varargin)
%SOLVE_RBSYSTEM solves the RB problem for a given parameter value
%
%   [UN] = SOLVE_RBSYSTEM(ROM, MU, N) given a parameter value MU and a 
%   dimension N (<= ROM.N), the function assembles the RB problem 
%   starting from the mu-independent arrays stored in ROM, and solves it. 
%   UN is the RB solution, given by the vector of coefficients associated 
%   to the reduced basis functions. 
%
%   [UN, UNH] = SOLVE_RBSYSTEM(ROM, MU, N) returns the high-fidelity 
%   representation UNH of UN.
%
%   [UN, UNH, BETAN] = SOLVE_RBSYSTEM(ROM, MU, N) returns the stability
%   factor BETAN of the RB problem

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 3 || isempty(N)
    N = ROM.N;
end

%% Evaluate Theta functions
[ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, ROM, varargin{:} );

INDX = 1:N;

%% Assemble RB problem
switch ROM.method
    case 'Galerkin'
        
        AN = 0 * ROM.ANq{1}(INDX,INDX);
        FN = 0 * ROM.FNq{1}(INDX);
        
        for q = 1 : ROM.Qa
            AN = AN + theta_a(q) * ROM.ANq{q}(INDX,INDX);
        end
        
        for q  = 1 : ROM.Qf
            FN = FN + theta_f(q) * ROM.FNq{q}(INDX);
        end
        
    case 'LeastSquares'
        
        AN = 0 * ROM.ANq{1,1}(INDX,INDX);
        FN = 0 * ROM.FNq{1,1}(INDX);
        
        for q1  = 1 : ROM.Qa
            for q2  = 1 : ROM.Qa
                AN = AN + theta_a(q1) * theta_a(q2) * ROM.ANq{q1,q2}(INDX,INDX);
            end
            
            for q2  = 1 : ROM.Qf
                FN = FN + theta_a(q1) * theta_f(q2) * ROM.FNq{q1,q2}(INDX);
            end
        end
        
end

%% Solve RB problem
uN   = AN \ FN;

%% Compute High-Fidelity representation of UN
if nargout > 1
    uNh                            = zeros(ROM.MESH.numNodes,1);
    uNh(ROM.MESH.internal_dof)     = ROM.V(:,INDX)*uN;
    
    switch ROM.model
        case 'ADR'
            uNh(ROM.MESH.Dirichlet_dof)    = ROM.u_D(ROM.MESH.nodes(:,ROM.MESH.Dirichlet_dof), mu);
            
        case 'CSM'
            for k = 1 : ROM.MESH.dim
                uNh(ROM.MESH.Dirichlet_dof_c{k}+(k-1)*ROM.MESH.numNodes)    = ROM.u_D{k}(ROM.MESH.nodes(:,ROM.MESH.Dirichlet_dof_c{k}), mu);
            end
    end
end

%% Compute stability factor of the RB problem
if nargout >= 3
    
    OPTS.tol     = 1e-5;
    OPTS.issym   = 1;
    OPTS.isreal  = 1;
    OPTS.disp    = 0;
    OPTS.maxit   = 2500;
    
    XN  = full(ROM.Xnorm(INDX, INDX));
    AN  = full(AN);
    
    switch ROM.method
        case 'Galerkin'
            betaN =  min(svd(XN^(-0.5)*(AN*XN^(-0.5))));
            
        case 'LeastSquares'
            betaN =  sqrt(eigs(AN,XN,1,'sm',OPTS));
    end
    
end


end