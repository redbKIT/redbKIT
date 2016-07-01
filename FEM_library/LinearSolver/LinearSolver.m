%LinearSolver solves linear systems Ax=b using direct and iterative solvers
% LinearSolver methods:
%    LinearSolver       - constructor
%    SetPreconditioner  - set preconditioner object
%    Solve              - solve linear system Ax = b
%
% LinearSolver properties:
%    M_options          - struct containing linear solver options
%    M_type             - linear solver type (automatically set by constructor)
%                         see M_options.type
%    M_precon           - Preconditioner object (only used in combination
%                         with iterative methods)

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

%{
properties
    %M_precon - see also Preconditioner.m
    M_precon;

    %M_options -
    %    type (mandatory): 'backslash', 'MUMPS', 'gmres', 'matlab_lu'
    %    mumps_reordering (only for type = 'MUMPS'):
    %               0 - Approximate Minimum Degree is used
    %               3 - SCOTCH (if available)
    %               4 - PORD (if available)
    %               5 - METIS (if available)
    %               7 - Automatic choice by MUMPS
    %    tol (only for type = 'gmres'): iterative solver tolerance
    %    maxit (only for type = 'gmres'): max number of iterations
    %    gmres_verbosity (only for type = 'gmres'): print convergence
    %               history each gmres_verbosity iterations
    M_options;

end
%}

classdef LinearSolver < handle
    
    properties (GetAccess = public, SetAccess = protected)
        M_type;
        M_options;
        M_precon;
        M_verbose;
        M_haveFactorization;
        M_solveTime;
    end
    
    properties (Access = private)
        M_L;
        M_U;
        M_perm;
        M_invperm;
    end
    
    methods
        
        %% Constructor
        function obj = LinearSolver( Options )
            
            obj.M_options   = Options;
            obj.M_type      = Options.type;
            obj.M_verbose   = false;
            obj.M_solveTime = 0;
            obj.M_haveFactorization = false;
            
        end
        
        %% Option Parser
        function obj = OptionParser( obj )
            % to be coded
        end
        
        %% SetPreconditioner
        function obj = SetPreconditioner( obj, Precon )
            %SetPreconditioner method
            %   LinearSolver.SetPreconditioner( Precon )
            
            obj.M_precon = Precon;
        end
        
        %% GetSolveTime
        function t = GetSolveTime( obj )
            t = obj.M_solveTime;
        end
        
        %% Solve
        function x = Solve( obj, A, b, x0 )
            %Solve method
            %   x = LinearSolver.Solve( A, b, x0 )
            
            if nargin < 4 || isempty(x0)
                x0 = zeros(length(b),1);
            end
            
            if obj.M_verbose
                fprintf('\n         Solving Linear System ...\n')
            end
            
            switch obj.M_type
                
                case 'backslash'
                    time_solve = tic;
                    
                    x = A \ b;
                    obj.M_solveTime = toc(time_solve);
                    
                case 'MUMPS'
                    time_solve = tic;
                    
                    % initialization of a matlab MUMPS structure
                    id     = initmumps;
                    id.SYM = 0;
                    % here JOB = -1, the call to MUMPS will initialize C
                    % and fortran MUMPS structure
                    id = dmumps(id);
                    % JOB = 6 means analysis + factorization + solve
                    id.JOB = 6;
                    
                    id.ICNTL(1:4) = -1; % no output
                    id.ICNTL(7)   = obj.M_options.mumps_reordering; % Typer of reordering
                    % set RHS
                    id.RHS = b;
                    
                    % Call Mumps
                    [id] = dmumps(id,A);
                    x = id.SOL;
                    
                    id.JOB = -2;
                    id = dmumps(id);
                    obj.M_solveTime = toc(time_solve);
                    
                case 'matlab_lu'
                    time_solve = tic;
                    
                    if  ~obj.M_haveFactorization
                        [obj.M_L , obj.M_U , obj.M_perm , q ]  = lu(A, 'vector');
                        obj.M_invperm             = 0*q ;
                        obj.M_invperm(q)          = 1:length(q);
                        obj.M_haveFactorization = true;
                    end
                    
                    x = obj.M_L \ b(obj.M_perm);
                    x = obj.M_U \ x;
                    x = x(obj.M_invperm);
                    
                    obj.M_solveTime = toc(time_solve);
                    
                case 'gmres'
                    
                    time_solve = tic;
                    
                    [x,flagITER,~,~,resvec] = my_gmres(A, b, 100, obj.M_options.tol, obj.M_options.maxit,...
                        @(r)obj.M_precon.Apply(r), [], x0, [1 obj.M_options.gmres_verbosity]);
                    
                    obj.M_solveTime = toc(time_solve);
                    
                    obj.M_precon.Clean();
                    
                    if flagITER == 0
                        fprintf('\nGmres converged in %d iterations\n',length(resvec));
                    else
                        fprintf('\n***Problems with the linear solver***\n');
                    end
                    
            end
            
            if obj.M_verbose
                fprintf('    Done in: %2.2e s  ', obj.M_solveTime);
            end
            
        end
        
    end
    
end