classdef LinearSolver < handle

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

    properties (GetAccess = public, SetAccess = protected)
        M_type;
        M_options;
        M_precon;
        M_verbose;
    end
   
    methods
        
        %% Constructor
        function obj = LinearSolver( Options )
            
            obj.M_options = Options;
            obj.M_type    = Options.type;
            obj.M_verbose = true;
            
        end
        
        %%
        function obj = SetPreconditioner( obj, Precon )
            obj.M_precon = Precon;
        end
        
        %%
        function x = Solve( obj, A, b, x0 )
            
            if nargin < 4 || isempty(x0)
                x0 = zeros(length(b),1);
            end

            %if obj.M_verbose
            %    fprintf('\n         Solving Linear System ...\n')
            %end

            time_solve = tic;
            
            switch obj.M_type
                
                case 'backslash'
                    
                    x = A \ b;
                    
                case 'MUMPS'
                    
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
                    time_solve = tic;
                    [id] = dmumps(id,A);
                    x = id.SOL;
                    time_solve = toc(time_solve);
                    
                    id.JOB = -2;
                    id = dmumps(id);
                    
                case 'gmres'
   
                    [x,flagITER,~,~,resvec] = my_gmres(A, b, [], obj.M_options.tol, obj.M_options.maxit,...
                        @(r)obj.M_precon.Apply(r), [], x0, [obj.M_verbose obj.M_options.gmres_verbosity]);
                    
            end
            
            time_solve = toc(time_solve);
            %if obj.M_verbose
            %    fprintf('    Done in: %2.2e s  ', time_solve);
            %end
            
        end
        
    end
    
end