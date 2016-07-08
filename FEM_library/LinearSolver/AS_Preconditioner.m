classdef AS_Preconditioner < Preconditioner & handle
%AS_PRECONDITIONER One-level additive schwarz preconditioner    
    
%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

    properties (GetAccess = public, SetAccess = protected)
        M_Restrictions;
        M_L;
        M_U;
        M_perm;
        M_invperm;
        M_mumps_ID;
        M_A_DD;
        M_twoLevel;
        M_A_coarse;
        M_R_coarse;
    end
    
    methods
        
        %% Constructor
        function obj = AS_Preconditioner( varargin )
            
            obj@Preconditioner( varargin{:} );
            
            if ~isfield(obj.M_options, 'coarse_level')
                obj.M_options.coarse_level = 'None';
            end
            
            % Parse Coarse Solver options
            switch obj.M_options.coarse_level
                case 'None'
                    obj.M_twoLevel = false;
                    
                case {'Aggregation', 'SmoothedAggregation'}
                    obj.M_twoLevel = true;
                    
                    if ~isfield(obj.M_options, 'coarse_num_aggregates')
                        obj.M_options.coarse_num_aggregates =  obj.M_options.num_subdomains;
                    end
                    
                    if strcmp(obj.M_options.coarse_level, 'SmoothedAggregation')
                        
                        if ~isfield(obj.M_options, 'coarse_smoother_iter')
                            obj.M_options.coarse_smoother_iter = 1;
                        end
                        
                        if ~isfield(obj.M_options, 'coarse_smoother_dumping')
                            obj.M_options.coarse_smoother_dumping = 1;
                        end
                    end
                    
                otherwise
                    error('AS_Preconditioner: invalid coarse_level type.');
                    
            end
            
        end
        
        %% Set Restriction Operators
        function obj = SetRestrictions(obj, Restriction_operators )
            obj.M_Restrictions = Restriction_operators;
        end
        
        %% Build preconditioner
        function obj = Build(obj, A )
            
            time_build = tic;
            
            switch obj.M_options.local_solver
                
                case 'matlab_lu'
                    
                    if ~obj.M_reuse || (obj.M_reuse && ~obj.M_isBuilt)
                        
                        A_DD = Composite(obj.M_options.num_subdomains);
                        for ii = 1:numel(A_DD)
                            A_DD{ii}             = A(obj.M_Restrictions{ii},obj.M_Restrictions{ii});
                        end

                        spmd(0,numel(A_DD))
                            [M_L , M_U , M_perm , q ]  = lu(A_DD, 'vector');
                            M_invperm          = 0*q ;
                            M_invperm(q)       = 1:length(q);
                        end
                        obj.M_L = M_L;
                        obj.M_U = M_U;
                        obj.M_perm    = M_perm;
                        obj.M_invperm = M_invperm;
                        
                        % coarse level
                        if obj.M_twoLevel
                            
                            if strcmp(obj.M_options.coarse_level, 'SmoothedAggregation')
                                D       = spdiags( spdiags(A, 0), 0, size(A,1), size(A,1));
                                R_coarse = obj.M_Restrictions{end}';
                                
                                for i = 1 : obj.M_options.coarse_smoother_dumping
                                    R_coarse = R_coarse - obj.M_options.coarse_smoother_dumping * ( D \ ( A*R_coarse ) );
                                end
                                
                                R_coarse       = R_coarse';
                            else
                                R_coarse = obj.M_Restrictions{end};
                            end
                            
                            obj.M_A_coarse = R_coarse * (A * R_coarse');
                            obj.M_R_coarse = R_coarse;
                            A_C = obj.M_A_coarse;
                            save A_C A_C;
                        end
                        
                        obj.M_isBuilt = true;
                        
                    end
                    
                case 'MUMPS'
                    
                    if ~obj.M_reuse || (obj.M_reuse && ~obj.M_isBuilt)
                        
                        A_DD = Composite(obj.M_options.num_subdomains);
                        for ii = 1:numel(A_DD)
                            A_DD{ii}   = A(obj.M_Restrictions{ii},obj.M_Restrictions{ii});
                        end
                        
                        reordering_alg = obj.M_options.mumps_reordering;
                        spmd(0,numel(A_DD))
                            mumps_ID            = initmumps;
                            mumps_ID.SYM        = 0;
                            mumps_ID            = dmumps(mumps_ID);
                            mumps_ID.JOB        = 4;
                            mumps_ID.ICNTL(1:4) = -1; % no output
                            mumps_ID.ICNTL(7)   = reordering_alg;
                            mumps_ID            = dmumps(mumps_ID, A_DD);
                        end
                        obj.M_A_DD     = A_DD;
                        obj.M_mumps_ID = mumps_ID;
                        clear mumps_ID;

                        % coarse level
                        if obj.M_twoLevel
                            if strcmp(obj.M_options.coarse_level, 'SmoothedAggregation')
                                D       = spdiags( spdiags(A, 0), 0, size(A,1), size(A,1));
                                R_coarse = obj.M_Restrictions{end}';
                                
                                for i = 1 : obj.M_options.coarse_smoother_dumping
                                    R_coarse = R_coarse - obj.M_options.coarse_smoother_dumping * ( D \ ( A*R_coarse ) );
                                end
                                
                                R_coarse       = R_coarse';
                            else
                                R_coarse = obj.M_Restrictions{end};
                            end
                            
                            obj.M_A_coarse = R_coarse * (A * R_coarse');
                            obj.M_R_coarse = R_coarse;
                            
                        end
                            
                        obj.M_isBuilt  = true;
                    end
                    
            end
            obj.M_BuildTime = toc(time_build);
            
        end
        
        %% Apply preconditioner
        function z = Apply(obj, r)

            switch obj.M_options.local_solver
                
                case 'matlab_lu'
                    
                    L =  obj.M_L;
                    U =  obj.M_U;
                    perm = obj.M_perm;
                    invperm = obj.M_invperm;
                    Restrictions = obj.M_Restrictions;
                    
                    % coarse level
                    if obj.M_twoLevel
                        warning ('off','all');
                        z  = obj.M_R_coarse'*(obj.M_A_coarse\(obj.M_R_coarse*r));
                        warning ('on','all');
                    else
                        z     = 0 *r;
                    end
                    
                    spmd(0,obj.M_options.num_subdomains)
                        zi = r(Restrictions{labindex});
                        
                        warning ('off','all');
                        zi_l = L \ zi(perm);
                        zi_l = U \ zi_l;
                        warning ('on','all');
                        
                        zi = 0*r;
                        zi(Restrictions{labindex}) = zi_l(invperm);
                    end
                    
                    for i = 1 : obj.M_options.num_subdomains
                        z = z + zi{i};
                    end
                    
                case 'MUMPS'
                    A_DD     = obj.M_A_DD;
                    mumps_ID = obj.M_mumps_ID;
                    Restrictions = obj.M_Restrictions;
                    
                    % coarse level
                    if obj.M_twoLevel
                        warning ('off','all');
                        z  = obj.M_R_coarse'*(obj.M_A_coarse\(obj.M_R_coarse*r));
                        warning ('on','all');
                    else
                        z     = 0 *r;
                    end
                    
                    spmd(0,obj.M_options.num_subdomains)
                        zi = 0*r;
                        
                        mumps_ID.JOB    = 3;
                        mumps_ID.RHS    = r(Restrictions{labindex});
                        mumps_ID        = dmumps(mumps_ID, A_DD);
                        
                        zi(Restrictions{labindex}) = mumps_ID.SOL;
                    end
                    
                    for i = 1 : obj.M_options.num_subdomains
                        z = z + zi{i};
                    end
                    
                    clear mumps_ID;
                    
            end
            
        end
        
        %% Clean preconditioner
        function obj = Clean( obj )
            
            switch obj.M_options.local_solver
                
                case 'MUMPS'
                    mumps_ID = obj.M_mumps_ID;
                    num_local = length(obj.M_Restrictions);
                    spmd(0,num_local)
                        mumps_ID.JOB    = -2;
                        mumps_ID        = dmumps(mumps_ID);
                    end
                    
                otherwise
                    % do nothing
            end
            
        end
        
    end
        
end