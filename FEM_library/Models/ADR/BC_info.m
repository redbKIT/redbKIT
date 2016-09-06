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

if ~isfield( MESH, 'dim' )
    MESH.dim = 2;
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
            type_Robin     = DATA.flag_robin{d};
            type_clamp_points = DATA.flag_clamp_points{d};
            
            if isempty(type_Dirichlet) && isempty(type_Neumann) && isempty(type_Pressure) && isempty(type_Robin) && isempty(type_clamp_points) 
                error(['No boundary conditions are imposed on component ', num2str(d)]);
            end
            
            %% Find Dirichlet dofs (if any)
            if ~isempty(type_Dirichlet)
                
                % Computes the Dirichlet dof of the domain
                nDir                    = length(type_Dirichlet);
                Dirichlet_side          = [];
                for kk = 1 : nDir
                    this_Dirichlet_side     = find(MESH.boundaries(bc_flag_row,:) == type_Dirichlet(kk));
                    this_Dirichlet_dof      = MESH.boundaries(1:MESH.numBoundaryDof, unique( this_Dirichlet_side ) );
                    MESH.DiriDof_CompFlag{d,kk}  = unique(this_Dirichlet_dof(:));
                    Dirichlet_side          = [Dirichlet_side, this_Dirichlet_side];
                end
                Dirichlet_side             = unique(Dirichlet_side);
                Dirichlet_dof              = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
                Dirichlet_dof              = unique( Dirichlet_dof(:) );
            else
                Dirichlet_dof = [];
            end
                
            dir_ringDofs = type_clamp_points;
            MESH.clamp_points{d} = dir_ringDofs;
                            
            if ~isempty(type_Dirichlet) || ~isempty(dir_ringDofs)  
                
                MESH.Dirichlet_dof_c{d}    = unique([Dirichlet_dof; dir_ringDofs]);
                MESH.internal_dof_c{d}     = setdiff([1:MESH.numNodes]',MESH.Dirichlet_dof_c{d});
                
            else
                MESH.internal_dof_c{d}  = [1:MESH.numNodes]';
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
            
            %% Find Robin boundaries (if any)
            if ~isempty(type_Robin)
                % Computes the Neumann dof of the domain
                nRob         = length(type_Robin);
                Robin_side = [];
                for k = 1 : nRob
                    Robin_side = [Robin_side,find(MESH.boundaries(bc_flag_row,:) == type_Robin(k))];
                end
                MESH.Robin_side{d} = unique(Robin_side);
            else
                MESH.Robin_side{d} = [];
            end
            
            
            %% Find Pressure boundaries (if any)
            if ~isempty(type_Pressure)
                % Computes the Pressure dof of the domain
                nPrex       = length(type_Pressure);
                Pressure_side = [];
                for kk = 1 : nPrex
                    this_Pressure_side = find(MESH.boundaries(bc_flag_row,:) == type_Pressure(kk));
                    MESH.Pressure_side_CompFlag{d,kk} = unique(this_Pressure_side);
                    Pressure_side = [Pressure_side, this_Pressure_side];
                end
                MESH.Pressure_side{d} = unique(Pressure_side);
            else
                MESH.Pressure_side{d} = [];
            end
        end
        
        %% Find DirichletNormal dofs (if any)
        type_DirichletNormal = DATA.flag_dirichletNormal;
        if ~isempty(type_DirichletNormal)
            
            i_C = [];
            j_C = [];
            v_C = [];
            
            for l = 1 : length(type_DirichletNormal)
                
                Dirichlet_side       = find(MESH.boundaries(bc_flag_row,:) == type_DirichletNormal(l));                
                Dirichlet_side       = unique(Dirichlet_side);
                Dirichlet_dof        = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
                MESH.DirichletNormal_dof{l}    = unique(Dirichlet_dof(:));
                MESH.DirichletNormal_N{l}      = mean( MESH.Normal_Faces(:, Dirichlet_side) , 2) ;
               
                one_N = ones(length(MESH.DirichletNormal_dof{l}), 1);
                
                if MESH.DirichletNormal_N{l}(1) == 0
                    
                    if MESH.DirichletNormal_N{l}(2) ~= 0
                        
                        MESH.DirichletNormal_dof{l} = MESH.DirichletNormal_dof{l} + MESH.numNpdes;
                        
                        i_C = [i_C; MESH.DirichletNormal_dof{l}];
                        j_C = [j_C; MESH.DirichletNormal_dof{l}-MESH.numNodes];
                        v_C = [v_C; -MESH.DirichletNormal_N{l}(1)/MESH.DirichletNormal_N{l}(2)*one_N];
                            
                        if MESH.dim == 3
                            
                            i_C = [i_C; MESH.DirichletNormal_dof{l}];
                            j_C = [j_C; MESH.DirichletNormal_dof{l}+MESH.numNodes];
                            v_C = [v_C; -MESH.DirichletNormal_N{l}(3)/MESH.DirichletNormal_N{l}(2)*one_N];
                            
                        end
                                                
                    elseif MESH.DirichletNormal_N{l}(3) ~= 0
                        
                        MESH.DirichletNormal_dof{l} = MESH.DirichletNormal_dof{l} + 2*MESH.numNpdes;
                        
                        for k = 1 : 2
                            i_C = [i_C; MESH.DirichletNormal_dof{l}];
                            j_C = [j_C; MESH.DirichletNormal_dof{l}-(3-k)*MESH.numNodes];
                            v_C = [v_C; -MESH.DirichletNormal_N{l}(k)/MESH.DirichletNormal_N{l}(3)*one_N];
                        end
                        
                    end
                    
                else
                    
                    for k = 2 : MESH.dim
                        i_C = [i_C; MESH.DirichletNormal_dof{l}];
                        j_C = [j_C; MESH.DirichletNormal_dof{l}+(k-1)*MESH.numNodes];
                        v_C = [v_C; -MESH.DirichletNormal_N{l}(k)/MESH.DirichletNormal_N{l}(1)*one_N];
                    end

                end
                
                
                MESH.internal_dof  = setdiff(MESH.internal_dof, MESH.DirichletNormal_dof{l});
            end
        
            i_C = [i_C; MESH.internal_dof];    
            j_C = [j_C; MESH.internal_dof];    
            v_C = [v_C; ones(length(MESH.internal_dof),1)];
            
            i_C = [i_C; MESH.Dirichlet_dof];    
            j_C = [j_C; MESH.Dirichlet_dof];    
            v_C = [v_C; ones(length(MESH.Dirichlet_dof),1)];
            
            MESH.DirichletNormal_R = sparse(i_C, j_C, v_C, MESH.dim*MESH.numNodes, MESH.dim*MESH.numNodes);
            
        else
            
            MESH.DirichletNormal_dof = [];
            MESH.DirichletNormal_N   = [];
            MESH.DirichletNormal_R   = speye( MESH.dim*MESH.numNodes, MESH.dim*MESH.numNodes );
            
        end        
        
        
    case 'CFD'
        
        MESH.Dirichlet_dof = [];
        MESH.internal_dof  = [];
        
        for d = 1 : MESH.dim
            
            type_Dirichlet = DATA.flag_dirichlet{d};
            type_Neumann   = DATA.flag_neumann{d};
            type_Pressure  = DATA.flag_pressure{d};
            type_rings     = DATA.flag_ring{d};
            
            if isempty(type_Dirichlet) && isempty(type_Neumann) && isempty(type_Pressure)
                error(['No boundary conditions are imposed on component ', num2str(d)]);
            end
            
            %% Find Dirichlet dofs (if any)
            if ~isempty(type_Dirichlet)
                
                % Computes the Dirichlet dof of the domain
                nDir                    = length(type_Dirichlet);
                Dirichlet_side          = [];
                for kk = 1 : nDir
                    this_Dirichlet_side     = find(MESH.boundaries(bc_flag_row,:) == type_Dirichlet(kk));
                    this_Dirichlet_dof      = MESH.boundaries(1:MESH.numBoundaryDof, unique( this_Dirichlet_side ) );
                    MESH.DiriDof_CompFlag{d,kk}  = unique(this_Dirichlet_dof(:));
                    Dirichlet_side          = [Dirichlet_side, this_Dirichlet_side];
                end
                Dirichlet_side             = unique(Dirichlet_side);
                Dirichlet_dof              = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
                Dirichlet_dof              = unique( Dirichlet_dof(:) );
            else
                Dirichlet_dof = [];
            end
                
            nRings = length(type_rings);
            dir_ringDofs = [];
            for j = 1 : nRings
                index        = find(MESH.rings(bc_flag_row,:) == type_rings(j));
                tmp          = MESH.rings(1:MESH.numRingsDof, index);
                dir_ringDofs = [dir_ringDofs unique(tmp(:))];
            end
            MESH.ringDofs{d} = dir_ringDofs;
                            
            if ~isempty(type_Dirichlet) || ~isempty(dir_ringDofs)  
                
                MESH.Dirichlet_dof_c{d}    = unique([Dirichlet_dof; dir_ringDofs]);
                MESH.internal_dof_c{d}     = setdiff([1:MESH.numNodes]',MESH.Dirichlet_dof_c{d});
                
            else
                MESH.internal_dof_c{d}  = [1:MESH.numNodes]';
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
                    this_Pressure_side = find(MESH.boundaries(bc_flag_row,:) == type_Pressure(kk));
                    MESH.Pressure_side_CompFlag{d,kk} = unique(this_Pressure_side);
                    Pressure_side = [Pressure_side, this_Pressure_side];
                end
                MESH.Pressure_side{d} = unique(Pressure_side);
            else
                MESH.Pressure_side{d} = [];
            end
            
        end
        
        %% Find Resistance boundaries (if any)
        type_resistance = DATA.flag_resistance;
        if ~isempty(type_resistance)
            % Find the Resistance dofs of the domain
            nPrex       = length(type_resistance);
            Resistance_side = [];
            for kk = 1 : nPrex
                this_Resistance_side = find(MESH.boundaries(bc_flag_row,:) == type_resistance(kk));
                MESH.Resistance_side_Flag{kk} = unique(this_Resistance_side);
                Resistance_side = [Resistance_side, this_Resistance_side];
            end
            MESH.Resistance_side = unique(Resistance_side);
        else
            MESH.Resistance_side = [];
        end
        
        %% Find Absorbing boundaries (if any)
        type_absorbing = DATA.flag_absorbing;
        if ~isempty(type_absorbing)
            % Find the Resistance dofs of the domain
            nPrex       = length(type_absorbing);
            Absorbing_side = [];
            for kk = 1 : nPrex
                this_Absorbing_side = find(MESH.boundaries(bc_flag_row,:) == type_absorbing(kk));
                MESH.Absorbing_side_Flag{kk} = unique(this_Absorbing_side);
                Absorbing_side = [Absorbing_side, this_Absorbing_side];
            end
            MESH.Absorbing_side = unique(Absorbing_side);
        else
            MESH.Absorbing_side = [];
        end

        MESH.internal_dof  = [MESH.internal_dof; MESH.dim*MESH.numNodes+[1:MESH.numVertices]' ];
        
    otherwise
        error('BC_info: unknown model')
        
end

MESH.bc_flag_row = bc_flag_row;

end

