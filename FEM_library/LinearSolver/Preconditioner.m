classdef Preconditioner < handle
%PRECONDITIONER abstract class for preconditioner object

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

    properties (GetAccess = public, SetAccess = protected)
        M_options;
        M_reuse;
        M_isBuilt;
        M_type;
        M_BuildTime;
    end
    
    methods
        
        %% Constructor
        function obj = Preconditioner(DATA)
                        
            obj.M_options   = DATA.Preconditioner;
            obj.M_isBuilt   = false;
            
            if ~isfield(obj.M_options, 'reuse')
                obj.M_reuse     = false;
            else
                obj.M_reuse     = obj.M_reuse;
            end      
            
            obj.M_type = obj.M_options.type;
            obj.M_BuildTime = 0;
            
        end
        
        %% GetBuildTime
        function t = GetBuildTime( obj )
            t = obj.M_BuildTime;
        end
        
        %% Build preconditioner
        function obj = Build(obj, A )
                       
        end
        
        %% Apply preconditioner
        function z = Apply(obj, r)
            z = r;
        end
        
        %% Clean preconditioner
        function obj = Clean( obj )
                       
        end
 
    end
        
end