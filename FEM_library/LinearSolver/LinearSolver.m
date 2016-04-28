classdef LinearSolver < handle
    
    %   This file is part of redbKIT.
    %   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
    %   Author: Federico Negri <federico.negri at epfl.ch>
    
    properties (GetAccess = public, SetAccess = protected)
        M_type;
        M_options;
        M_precon;
        M_verbose;
        M_solveTime;
    end
    
    methods
        
        %% Constructor
        function obj = LinearSolver( Options )
            
            obj.M_options   = Options;
            obj.M_type      = Options.type;
            obj.M_verbose   = false;
            obj.M_solveTime = 0;
            
        end
        
        %% Option Parser
        function obj = OptionParser( obj )
           % to be coded 
        end
        
        %% SetPreconditioner
        function obj = SetPreconditioner( obj, Precon )
            obj.M_precon = Precon;
        end
        
        %% GetSolveTime
        function t = GetSolveTime( obj )
            t = obj.M_solveTime;
        end
        
        %% Solve
        function x = Solve( obj, A, b, x0 )
            
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
                    
                case 'gmres'
                    
                    time_solve = tic;
                    
                    [x,flagITER,~,~,resvec] = my_gmres(A, b, 100, obj.M_options.tol, obj.M_options.maxit,...
                        @(r)obj.M_precon.Apply(r), [], x0, [1 obj.M_options.gmres_verbosity]);
                    
                    obj.M_solveTime = toc(time_solve);
                    
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