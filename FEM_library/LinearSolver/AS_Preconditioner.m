classdef AS_Preconditioner < Preconditioner & handle
    
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
    end
    
    methods
        
        %% Constructor
        function obj = AS_Preconditioner( varargin )
            
            obj@Preconditioner( varargin{:} );
            
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
                        
                        A_DD = Composite(length(obj.M_Restrictions));
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
                        
                        obj.M_perm = M_perm;
                        obj.M_invperm = M_invperm;
                        obj.M_isBuilt = true;
                    end
                    
                case 'MUMPS'
                    
                    if ~obj.M_reuse || (obj.M_reuse && ~obj.M_isBuilt)
                        
                        A_DD = Composite(length(obj.M_Restrictions));
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
                     
                    spmd(0,length(Restrictions))
                        zi = r(Restrictions{labindex});
                        
                        warning ('off','all');
                        zi_l = L \ zi(perm);
                        zi_l = U \ zi_l;
                        warning ('on','all');
                        
                        zi = 0*r;
                        zi(Restrictions{labindex}) = zi_l(invperm);
                    end
                    
                    z     = 0 *r;
                    for i = 1 : length(obj.M_Restrictions)
                        z = z + zi{i};
                    end
                    
                case 'MUMPS'
                    A_DD     = obj.M_A_DD;
                    mumps_ID = obj.M_mumps_ID;
                    Restrictions = obj.M_Restrictions;
                    spmd(0,length(Restrictions))
                        zi = 0*r;
                        
                        obj.M_mumps_ID.JOB    = 3;
                        obj.M_mumps_ID.RHS    = r(Restrictions{labindex});
                        obj.M_mumps_ID        = dmumps(mumps_ID, A_DD);
                        
                        zi(Restrictions{labindex}) = mumps_ID.SOL;
                    end
                    
                    z     = 0 *r;
                    for i = 1 : length(obj.M_Restrictions)-1
                        z = z + zi{i};
                    end
                    
            end
            
        end
        
    end
        
end