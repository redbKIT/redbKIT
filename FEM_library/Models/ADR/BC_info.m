function [MESH] = BC_info(MESH, DATA, model)
%BC_INFO update MESH struct fields with Boundary conditions information
%
%   [MESH] = BC_INFO(MESH, DATA) given a MESH struct (see below for the
%   required fields) and a DATA structure containing the problem setup,
%   updates MESH adding all the fields required to impose boundary
%   conditions.
%
%   The input MESH should contain the following fields:
%   nodes, elements, boundaries, numNodes
%
%   The output MESH contains the following additional fields:
%   internal_dofs, Dirichlet_dofs, Neumann_side, Robin_side

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 3 || isempty(model)
    model = 'ADR';
end

%% Parse MESH struct to check if the the required fields are available
MESH_parser = inputParser;
MESH_parser.KeepUnmatched = true;

if verLessThan('matlab', '8.2')
    addParamValue(MESH_parser,'nodes',0);
    addParamValue(MESH_parser,'boundaries',0);
    addParamValue(MESH_parser,'elements',0);
    addParamValue(MESH_parser,'numNodes',0);
    addParamValue(MESH_parser,'dim',0);
else
    addParameter(MESH_parser,'nodes',0);
    addParameter(MESH_parser,'boundaries',0);
    addParameter(MESH_parser,'elements',0);
    addParameter(MESH_parser,'numNodes',0);
    addParameter(MESH_parser,'dim',0);
end

parse(MESH_parser,MESH);

if ~isempty(MESH_parser.UsingDefaults)
    msg = sprintf('MESH data structure does not contain the following fields: \n');
    for i = 1 : length(MESH_parser.UsingDefaults)
        msg = [msg, ' ', MESH_parser.UsingDefaults{i}];
    end
    error(msg)
end

switch MESH.dim
    case 2
        bc_flag_row = 5;
    case 3
        bc_flag_row = 12;
end

switch model
    case 'ADR'
        
        type_Dirichlet = DATA.flag_dirichlet;
        type_Neumann   = DATA.flag_neumann;
        type_Robin     = DATA.flag_robin;
        
        if isempty(type_Dirichlet) && isempty(type_Neumann) && isempty(type_Robin)
            error('No boundary conditions are imposed');
        end
        
        %% Find Dirichlet dofs (if any)
        if ~isempty(type_Dirichlet)
            
            % Computes the Dirichlet dof of the domain
            nDir                    = length(type_Dirichlet);
            Dirichlet_side          = [];
            flag_Dirichlet_vertices = [];
            for k = 1 : nDir
                Dirichlet_side          = [Dirichlet_side,find(MESH.boundaries(bc_flag_row,:) == type_Dirichlet(k))];
                flag_Dirichlet_vertices = [flag_Dirichlet_vertices,type_Dirichlet(k)];
            end
            Dirichlet_side             = unique(Dirichlet_side);
            Dirichlet_dof              = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
            MESH.Dirichlet_dof         = unique(Dirichlet_dof(:));
            MESH.internal_dof          = setdiff([1:MESH.numNodes],MESH.Dirichlet_dof);
            
        else
            MESH.internal_dof  = 1:MESH.numNodes;
            MESH.Dirichlet_dof = [];
        end
        
        
        %% Find Neumann boundaries (if any)
        if ~isempty(type_Neumann)
            
            % Computes the Neumann dof of the domain
            nNeu         = length(type_Neumann);
            Neumann_side = [];
            for k = 1 : nNeu
                Neumann_side = [Neumann_side,find(MESH.boundaries(bc_flag_row,:) == type_Neumann(k))];
            end
            MESH.Neumann_side = unique(Neumann_side);
            
        else
            MESH.Neumann_side = [];
        end
        
        %% Find Robin boundaries (if any)
        if ~isempty(type_Robin)
            
            % Computes the Robin dof of the domain
            nRob       = length(type_Robin);
            Robin_side = [];
            for k = 1 : nRob
                Robin_side = [Robin_side,find(MESH.boundaries(bc_flag_row,:) == type_Robin(k))];
            end
            MESH.Robin_side = unique(Robin_side);
            
        else
            MESH.Robin_side = [];
        end
        
    case 'CSM'
        
        MESH.Dirichlet_dof = [];
        MESH.internal_dof  = [];
            
        for d = 1 : MESH.dim
            
            type_Dirichlet = DATA.flag_dirichlet{d};
            type_Neumann   = DATA.flag_neumann{d};
            type_Pressure  = DATA.flag_pressure{d};
            
            if isempty(type_Dirichlet) && isempty(type_Neumann) && isempty(type_Pressure)
                error(['No boundary conditions are imposed on component ', num2str(d)]);
            end
            
            %% Find Dirichlet dofs (if any)
            if ~isempty(type_Dirichlet)
                
                % Computes the Dirichlet dof of the domain
                nDir                    = length(type_Dirichlet);
                Dirichlet_side          = [];
                flag_Dirichlet_vertices = [];
                for kk = 1 : nDir
                    Dirichlet_side          = [Dirichlet_side,find(MESH.boundaries(bc_flag_row,:) == type_Dirichlet(kk))];
                    flag_Dirichlet_vertices = [flag_Dirichlet_vertices,type_Dirichlet(kk)];
                end
                Dirichlet_side             = unique(Dirichlet_side);
                Dirichlet_dof              = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
                MESH.Dirichlet_dof_c{d}    = unique(Dirichlet_dof(:));
                MESH.internal_dof_c{d}     = setdiff([1:MESH.numNodes]',MESH.Dirichlet_dof_c{d});
                
            else
                MESH.internal_dof_c{d}  = 1:MESH.numNodes;
                MESH.Dirichlet_dof_c{d} = [];
            end
            
            
            MESH.Dirichlet_dof = [MESH.Dirichlet_dof;  (d-1)*MESH.numNodes+MESH.Dirichlet_dof_c{d}];
            MESH.internal_dof  = [MESH.internal_dof; (d-1)*MESH.numNodes+MESH.internal_dof_c{d}];
            
            %% Find Neumann boundaries (if any)
            if ~isempty(type_Neumann)
                % Computes the Neumann dof of the domain
                nNeu         = length(type_Neumann);
                Neumann_side = [];
                for k = 1 : nNeu
                    Neumann_side = [Neumann_side,find(MESH.boundaries(bc_flag_row,:) == type_Neumann(k))];
                end
                MESH.Neumann_side{d} = unique(Neumann_side);
            else
                MESH.Neumann_side{d} = [];
            end
            
            
            %% Find Pressure boundaries (if any)
            if ~isempty(type_Pressure)
                % Computes the Pressure dof of the domain
                nPrex       = length(type_Pressure);
                Pressure_side = [];
                for kk = 1 : nPrex
                    Pressure_side = [Pressure_side,find(MESH.boundaries(bc_flag_row,:) == type_Pressure(kk))];
                end
                MESH.Pressure_side{d} = unique(Pressure_side);
            else
                MESH.Pressure_side{d} = [];
            end
        end
        
        
    case 'CFD'
        
        MESH.Dirichlet_dof = [];
        MESH.internal_dof  = [];
        
        for d = 1 : MESH.dim
            
            type_Dirichlet = DATA.flag_dirichlet{d};
            type_Neumann   = DATA.flag_neumann{d};
            type_Pressure  = DATA.flag_pressure{d};
            
            if isempty(type_Dirichlet) && isempty(type_Neumann) && isempty(type_Pressure)
                error(['No boundary conditions are imposed on component ', num2str(d)]);
            end
            
            %% Find Dirichlet dofs (if any)
            if ~isempty(type_Dirichlet)
                
                % Computes the Dirichlet dof of the domain
                nDir                    = length(type_Dirichlet);
                Dirichlet_side          = [];
                flag_Dirichlet_vertices = [];
                for kk = 1 : nDir
                    Dirichlet_side          = [Dirichlet_side,find(MESH.boundaries(bc_flag_row,:) == type_Dirichlet(kk))];
                    flag_Dirichlet_vertices = [flag_Dirichlet_vertices,type_Dirichlet(kk)];
                end
                Dirichlet_side             = unique(Dirichlet_side);
                Dirichlet_dof              = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
                MESH.Dirichlet_dof_c{d}    = unique(Dirichlet_dof(:));
                MESH.internal_dof_c{d}     = setdiff([1:MESH.numNodes]',MESH.Dirichlet_dof_c{d});
                
            else
                MESH.internal_dof_c{d}  = 1:MESH.numNodes;
                MESH.Dirichlet_dof_c{d} = [];
            end
            
            
            MESH.Dirichlet_dof = [MESH.Dirichlet_dof;  (d-1)*MESH.numNodes+MESH.Dirichlet_dof_c{d}];
            MESH.internal_dof  = [MESH.internal_dof; (d-1)*MESH.numNodes+MESH.internal_dof_c{d}];
            
            %% Find Neumann boundaries (if any)
            if ~isempty(type_Neumann)
                % Computes the Neumann dof of the domain
                nNeu         = length(type_Neumann);
                Neumann_side = [];
                for k = 1 : nNeu
                    Neumann_side = [Neumann_side,find(MESH.boundaries(bc_flag_row,:) == type_Neumann(k))];
                end
                MESH.Neumann_side{d} = unique(Neumann_side);
            else
                MESH.Neumann_side{d} = [];
            end
            
            %% Find Pressure boundaries (if any)
            if ~isempty(type_Pressure)
                % Computes the Pressure dof of the domain
                nPrex       = length(type_Pressure);
                Pressure_side = [];
                for kk = 1 : nPrex
                    Pressure_side = [Pressure_side,find(MESH.boundaries(bc_flag_row,:) == type_Pressure(kk))];
                end
                MESH.Pressure_side{d} = unique(Pressure_side);
            else
                MESH.Pressure_side{d} = [];
            end
            
        end
        
        MESH.internal_dof  = [MESH.internal_dof; MESH.dim*MESH.numNodes+[1:MESH.numVertices]' ];
        
    otherwise
        error('BC_info: unknown model')
        
end

MESH.bc_flag_row = bc_flag_row;

end

