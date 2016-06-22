%NEWMARK_TimeAdvance class to Handle the time advancing scheme based on 
%Newmark scheme

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

classdef Newmark_TimeAdvance < handle
    
    properties (GetAccess = public, SetAccess = protected)
        M_U;
        M_dU;
        M_d2U;
        M_Csi;
        M_beta;
        M_gamma;
        M_timestep;
        M_stateSize;
    end
    
    methods (Access = public)
        
        %% Constructor
        function obj = Newmark_TimeAdvance( beta, gamma, timestep )
            obj.M_beta      = beta;
            obj.M_gamma     = gamma;
            obj.M_timestep  = timestep;
            obj.PrintCoefficients();
        end
        
        %% Initialize
        function obj = Initialize( obj, U, dU, d2U )
            if size(U,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of U')
            end
            obj.M_stateSize = size(U, 1);
            obj.M_U   = U;
            
            if size(dU,1) ~= obj.M_stateSize || size(dU,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of dU')
            end
            obj.M_dU  = dU;
            
            if size(d2U,1) ~= obj.M_stateSize || size(d2U,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of d2U')
            end
            
            obj.M_d2U = d2U;
            
            obj.M_Csi = 1 / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
        end
        
        %% Update
        function obj = Update( obj, U_np1 )
            
            d2U_np1        = 1 / ( obj.M_beta * obj.M_timestep^2) * U_np1 - obj.M_Csi;
            obj.M_dU       = obj.M_dU + obj.M_timestep * (obj.M_gamma * d2U_np1 + (1 - obj.M_gamma) * obj.M_d2U) ;
            obj.M_d2U      = d2U_np1;
            obj.M_U        = U_np1;
            
            obj.M_Csi = 1 / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
            
        end
        
        %% RhsContribute
        function Csi = RhsContribute( obj )
            
            Csi = obj.M_Csi;
            
        end
        
        %% Mass coefficient
        function c = MassCoefficient( obj )
            
            c = 1 / ( obj.M_beta * obj.M_timestep^2) ;
            
        end
        
        %% PrintCoefficients
        function obj = PrintCoefficients( obj )
            
            fprintf('\n=======================================');
            fprintf('\nNewmark Time Advance Class');
            fprintf('\n  beta     = %1.3f', obj.M_beta);
            fprintf('\n  gamma    = %2.3f', obj.M_gamma);
            fprintf('\n  timestep = %1.3e', obj.M_timestep);
            fprintf('\n=======================================\n');
            
        end
    end
    
end