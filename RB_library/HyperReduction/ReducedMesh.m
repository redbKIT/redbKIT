classdef ReducedMesh < handle
%REDUCEDMESH Class
    
%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

    properties (GetAccess = public, SetAccess = protected)
        M_MESH;
        M_Red_Mesh;
        M_fem;
        M_node_to_element;
        M_node_to_boundary;
        M_PDEtype;
        M_DoFsList;
        M_NodesList;
        M_ReducedElements;
        M_ReducedBoundaries;
    end
    
    properties (GetAccess = private, SetAccess = protected)
        M_numInternalDoFs;
        M_numDoFs;
    end
    
    methods
        
        %% Constructor
        function obj = ReducedMesh( MESH, fem, PDE_type )
            
            obj.M_MESH = MESH;
            obj.M_fem    = fem;
            
            if nargin < 3
                obj.M_PDEtype = 'ADR';
            else
                obj.M_PDEtype = PDE_type;
            end
            
            [ ~, obj.M_node_to_element, obj.M_node_to_boundary ] = compute_adjacency_elements(obj.M_MESH.nodes, ...
                obj.M_MESH.elements, obj.M_MESH.dim, obj.M_MESH.boundaries, obj.M_fem);
            
            obj.M_DoFsList = [];
            obj.M_NodesList        = [];
            
            obj.M_numInternalDoFs  =  length(obj.M_MESH.internal_dof); 
            
        end
        
        %% AppendInternalDoFs
        function obj = AppendInternalDoFs(obj, Dofs)
            
            %tmp         = zeros(obj.M_numDoFs,1);
            %tmp(obj.M_MESH.internal_dof(Dofs)) = 1;
            
            if size(obj.M_MESH.internal_dof(Dofs), 1) > 1
                obj.M_DoFsList = [obj.M_DoFsList obj.M_MESH.internal_dof(Dofs)'];
            else
                obj.M_DoFsList = [obj.M_DoFsList obj.M_MESH.internal_dof(Dofs)];
            end
            
        end
        
        %% AppendInternalDoFs_Vectorized
        function obj = AppendInternalDoFs_Vectorized(obj, Dofs)
            
            tmp        = sparse(obj.M_numInternalDoFs*obj.M_numInternalDoFs,1);
            tmp(Dofs)  = 1;
            
            tmp        = reshape(tmp, obj.M_numInternalDoFs, obj.M_numInternalDoFs);

            [row, col] = find(tmp);
            
            if size(obj.M_MESH.internal_dof, 1) > 1
                obj.M_DoFsList = [obj.M_DoFsList obj.M_MESH.internal_dof(unique([row; col]))'];
            else
                obj.M_DoFsList = [obj.M_DoFsList obj.M_MESH.internal_dof(unique([row; col]))];
            end
            
        end
        
        %% DoFs To Nodes
        function obj = DoFsToNodes( obj )
            
            switch obj.M_PDEtype
                
                case 'ADR'
                    
                    obj.M_NodesList = obj.M_DoFsList;
                    
                case 'CSM'
                    
                    obj.M_NodesList = obj.M_DoFsList - floor( ( obj.M_DoFsList - 1 ) ./ obj.M_MESH.numNodes ) * obj.M_MESH.numNodes;
                    
                case 'CFD'
                    
                    obj.M_NodesList = obj.M_DoFsList - floor( ( obj.M_DoFsList - 1 ) ./ obj.M_MESH.numNodes ) * obj.M_MESH.numNodes;          
            end
            
        end
        
        %% FindReducedElements
        function obj = FindReducedElements( obj )
            
            obj = obj.DoFsToNodes( );
            
            % Find Volume Elements
            Red_elem = [];
            for i = 1 : length(obj.M_NodesList)
                Red_elem  = [Red_elem obj.M_node_to_element{obj.M_NodesList(i)}];
            end
            obj.M_ReducedElements = unique(Red_elem);
            
            % Find Boundary Elements
            Attached_nodes      = obj.M_MESH.elements(1:obj.M_MESH.numElemDof,obj.M_ReducedElements);
            Attached_nodes      = unique(Attached_nodes(:));

            Red_boundary = [];
            for i = 1 : length(Attached_nodes)
                Red_boundary  = [Red_boundary obj.M_node_to_boundary{Attached_nodes(i)}];
            end
            obj.M_ReducedBoundaries = unique(Red_boundary);
                                    
        end
        
        %% Export Reduced Mesh to Vtk file
        function ExportToVtk( obj, fig_folder, filename )
            
            if nargin < 2 || isempty(fig_folder)
                fig_folder = '';
            else
                [~,~,~] = mkdir(fig_folder);
            end
            
            if nargin < 3 || isempty(filename)
                filename = '';
            end
            
            % Save reduced mesh to vtk for visualization
            ADR_export_solution(obj.M_MESH.dim, ones(obj.M_MESH.numVertices,1), ...
                obj.M_MESH.vertices, obj.M_MESH.elements, [fig_folder,filename,'_ReferenceMesh']);
            
            ADR_export_solution(obj.M_MESH.dim, ones(obj.M_MESH.numVertices,1), ...
                obj.M_MESH.vertices, obj.M_MESH.elements(:,obj.M_ReducedElements), [fig_folder,filename,'_ReducedMesh']);
            
        end
        
        %% Build ReducedMesh Data Structure
        function obj = Build( obj, DATA )
            
             obj = obj.FindReducedElements( );
             
            [ obj.M_Red_Mesh ] = buildMESH( obj.M_MESH.dim, obj.M_MESH.elements, ...
                obj.M_MESH.vertices, obj.M_MESH.boundaries, obj.M_fem, obj.M_MESH.quad_order, ...
                DATA, obj.M_PDEtype, [], obj.M_ReducedElements, obj.M_ReducedBoundaries );
            
        end
        
    end
        
end