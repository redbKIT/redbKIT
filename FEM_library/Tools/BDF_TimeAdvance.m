%BDF_TIMEADVANCE class to Handle the time advancing scheme based on BDF
%schemes

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

classdef BDF_TimeAdvance < handle
    
    properties (GetAccess = public, SetAccess = protected)
        M_currentOrder;
        M_order;
        M_states;
        M_stateSize;
        M_coeffBDF;
        M_coeffEXT;
    end
    
    methods (Access = public)
        
        %% Constructor
        function obj = BDF_TimeAdvance( order )
            obj.M_currentOrder = 1;
            obj.M_order        = order;
            obj.UpdateCoefficients();
            obj.PrintCoefficients();
        end
        
        %% Initialize
        function obj = Initialize( obj, state )
            if size(state,2) > 1
                error('BDF_TimeAdvance.Initialize: incorrectSize of state')
            end
            obj.M_stateSize = size(state, 1);
            obj.M_states{1} = state;
        end
        
        %% Append
        function obj = Append( obj, state )
            if size(state,2) > 1 || size(state,1)~= obj.M_stateSize
                error('BDF_TimeAdvance.Append: incorrectSize of state')
            end
            
            if obj.M_currentOrder == obj.M_order
                
                if obj.M_currentOrder > 1
                    
                    for i = 1 : obj.M_currentOrder-1
                        obj.M_states{i} = obj.M_states{i+1};
                    end
                    
                end
                obj.M_states{obj.M_currentOrder} = state;
                
            else
                obj.M_states{obj.M_currentOrder+1} = state;
                obj.M_currentOrder = obj.M_currentOrder + 1;
                obj.UpdateCoefficients();
                obj.PrintCoefficients();
            end
            
        end
        
        %% Extrapolate
        function u_ext = Extrapolate( obj )
            
            u_ext = zeros(obj.M_stateSize,1);
            for i = 1 : obj.M_currentOrder
                u_ext = u_ext + obj.M_states{obj.M_currentOrder+1-i} * obj.M_coeffEXT(i);
            end
            
        end
        
        %% Rhs Contribute
        function u_rhs = RhsContribute( obj )
            
            u_rhs = zeros(obj.M_stateSize,1);
            for i = 1 : obj.M_currentOrder
                u_rhs = u_rhs + obj.M_states{obj.M_currentOrder+1-i} * obj.M_coeffBDF(i+1);
            end
            
        end
        
        %% GetCoefficientDerivative
        function alpha = GetCoefficientDerivative( obj )
            
            alpha = obj.M_coeffBDF(1);
            
        end
        
        %% PrintCoefficients
        function obj = PrintCoefficients( obj )
            
            fprintf('\nBDF Time Advance Class: M_currentOrder is %d', obj.M_currentOrder);
            fprintf('\n  BDF Coefficients:');
            for i = 1 : length(obj.M_coeffBDF)
                fprintf(' %2.2f', obj.M_coeffBDF(i));
            end
            fprintf('\n  Extrapolation Coefficients:');
            for i = 1 : length(obj.M_coeffEXT)
                fprintf(' %2.2f', obj.M_coeffEXT(i));
            end
            fprintf('\n');
            
        end
        
    end
    
    
    methods (Access = private)
        %% UpdateCoefficients
        function obj = UpdateCoefficients( obj )
            
            switch obj.M_currentOrder
                
                case 1
                    obj.M_coeffBDF = [1; 1];
                    obj.M_coeffEXT = [1];
                    
                case 2
                    obj.M_coeffBDF = [3/2; 2; -1/2];
                    obj.M_coeffEXT = [2; -1];
                    
                case 3
                    obj.M_coeffBDF = [11/6; 3; -3/2; 1/3];
                    obj.M_coeffEXT = [3; -3; 1];
                    
                case 4
                    obj.M_coeffBDF = [25/12; 4; -3; 4/3; -1/4];
                    obj.M_coeffEXT = [4; -6; 4; -1];
                    
                otherwise
                    error('Unimplemented BDF scheme. Only first, second, 3th and 4th order schemes are available.');
            end
            
        end
        
    end
    
end