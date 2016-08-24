classdef AS_Preconditioner_Serial < Preconditioner & handle
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
        M_A;
    end
    
    methods
        
        %% Constructor
        function obj = AS_Preconditioner_Serial( varargin )
            
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
            obj.M_A = A;
            switch obj.M_options.local_solver
                
                case 'matlab_lu'
                    
                    if ~obj.M_reuse || (obj.M_reuse && ~obj.M_isBuilt)
                        
                        A_DD = cell(obj.M_options.num_subdomains,1);
                        for ii = 1 : length(A_DD) 
                            A_DD{ii}             = A(obj.M_Restrictions{ii},obj.M_Restrictions{ii});
                        end

                        for i = 1 : length(A_DD) 
                            [M_L{i} , M_U{i} , M_perm{i} , q ]  = lu(A_DD{i}, 'vector');
                            M_invperm{i}          = 0*q ;
                            M_invperm{i}(q)       = 1:length(q);
                        end
                        
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
                            
                            i = obj.M_options.num_subdomains + 1;
                            [M_L{i} , M_U{i} , M_perm{i} , q ]  = lu(obj.M_A_coarse, 'vector');
                            M_invperm{i}          = 0*q ;
                            M_invperm{i}(q)       = 1:length(q);
                        end
                        
                        obj.M_L = M_L;
                        obj.M_U = M_U;
                        obj.M_perm    = M_perm;
                        obj.M_invperm = M_invperm;
                        
                        obj.M_isBuilt = true;
                        
                    end
                    
                case 'MUMPS'
                    
                    if ~obj.M_reuse || (obj.M_reuse && ~obj.M_isBuilt)
                        
                        A_DD = cell(obj.M_options.num_subdomains,1);
                        for ii = 1 : length(A_DD) 
                            A_DD{ii}             = A(obj.M_Restrictions{ii},obj.M_Restrictions{ii});
                        end
                        
                        reordering_alg = obj.M_options.mumps_reordering;
                        mumps_ID       = cell(obj.M_options.num_subdomains,1);
                        for i = 1 : length(A_DD) 
                            mumps_ID{i}            = initmumps;
                            mumps_ID{i}.SYM        = 0;
                            mumps_ID{i}            = dmumps(mumps_ID{i});
                            mumps_ID{i}.JOB        = 4;
                            mumps_ID{i}.ICNTL(1:4) = -1; % no output
                            mumps_ID{i}.ICNTL(7)   = reordering_alg;
                            mumps_ID{i}            = dmumps(mumps_ID{i}, A_DD{i});
                        end
                        obj.M_A_DD     = A_DD;
                        obj.M_mumps_ID = mumps_ID;    
                        obj.M_isBuilt  = true;
                        clear mumps_ID;
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
                    
                    % coarse level (additive)
                    if obj.M_twoLevel
                        
                       i = obj.M_options.num_subdomains + 1;
                       
                       %warning ('off','all');
                       %z  = obj.M_R_coarse'*(obj.M_A_coarse\(obj.M_R_coarse*r));
                       %warning ('on','all');
                       
                       w = obj.M_R_coarse*r;
                       
                       warning ('off','all');
                       zi_l = L{i} \ w(perm{i});
                       zi_l = U{i} \ zi_l;
                       warning ('on','all');
                       
                       z = obj.M_R_coarse' * ( zi_l(invperm{i}) );                       
                    else
                       z     = 0 *r;
                    end
                    
                    zi = cell(obj.M_options.num_subdomains,1);
                    for i =  1 : obj.M_options.num_subdomains
                        zi{i} = r(Restrictions{i});
                        
                        warning ('off','all');
                        zi_l = L{i} \ zi{i}(perm{i});
                        zi_l = U{i} \ zi_l;
                        warning ('on','all');
                        
                        zi{i} = 0*r;
                        zi{i}(Restrictions{i}) = zi_l(invperm{i});
                    end
                    
                    %z     = 0 *r;
                    for i = 1 : obj.M_options.num_subdomains
                        z = z + zi{i};
                    end
                    
                    % coarse level (multiplicative)
                    %if obj.M_twoLevel
                    %    %warning ('off','all');
                    %    z  = z + obj.M_R_coarse'*(obj.M_A_coarse\(obj.M_R_coarse*(r-obj.M_A*z)));
                    %    %warning ('on','all');
                    %end
                    
                case 'MUMPS'
                    Restrictions = obj.M_Restrictions;
                    
                    % coarse level
                    if obj.M_twoLevel
                        warning ('off','all');
                        z  = obj.M_R_coarse'*(obj.M_A_coarse\(obj.M_R_coarse*r));
                        warning ('on','all');
                    else
                        z     = 0 *r;
                    end
                    
                    for i =  1 : obj.M_options.num_subdomains
                        zi{i} = 0*r;
                        
                        obj.M_mumps_ID{i}.JOB    = 3;
                        obj.M_mumps_ID{i}.RHS    = r(Restrictions{i});
                        obj.M_mumps_ID{i}        = dmumps(obj.M_mumps_ID{i}, obj.M_A_DD{i});
                        
                        zi{i}(Restrictions{i}) = obj.M_mumps_ID{i}.SOL;
                    end
                    
                    for i = 1 : obj.M_options.num_subdomains
                        z = z + zi{i};
                    end
                                        
            end
            
        end
        
        %% Clean preconditioner
        function obj = Clean( obj )
            
            switch obj.M_options.local_solver
                
                case 'MUMPS'
                    mumps_ID = obj.M_mumps_ID;
                    num_local = length(obj.M_Restrictions);
                    for i = 1 : length( mumps_ID )
                        mumps_ID{i}.JOB    = -2;
                        mumps_ID{i}        = dmumps(mumps_ID{i});
                    end
                    
                otherwise
                    % do nothing
            end
            
        end
        
    end
        
end