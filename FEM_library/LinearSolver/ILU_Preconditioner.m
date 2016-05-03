classdef ILU_Preconditioner < Preconditioner & handle
    
%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

    properties (GetAccess = public, SetAccess = protected)
        M_L;
        M_U;
        M_P;
        M_Setup;
    end
    
    methods
        
        %% Constructor
        function obj = ILU_Preconditioner( varargin )
            
            obj@Preconditioner( varargin{:} );
            
        end
  
        %% Build preconditioner
        function obj = Build(obj, A )
            
            if ~obj.M_reuse || (obj.M_reuse && ~obj.M_isBuilt)
                
                time_build = tic;
                obj.M_Setup.type    = obj.M_options.ILU_type;
                obj.M_Setup.droptol = obj.M_options.ILU_droptol;
                obj.M_Setup.milu    = 'row';
                
                [obj.M_L, obj.M_U, obj.M_P] = ilu(A, obj.M_Setup);
                
                obj.M_BuildTime = toc(time_build);
                obj.M_isBuilt = true;
                
            end
            
        end
        
        %% Apply preconditioner
        function z = Apply(obj, r)
    
            r = obj.M_P * r;
            
            z = obj.M_U \ ( obj.M_L \ r);
            
        end
        
    end
        
end