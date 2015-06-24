function [MESH] = BC_info(MESH, DATA)
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

%% Parse MESH struct to check if the the required fields are available
MESH_parser = inputParser;
MESH_parser.KeepUnmatched = true;

if verLessThan('matlab', '8.2')
      addParamValue(MESH_parser,'nodes',0);
      addParamValue(MESH_parser,'boundaries',0);
      addParamValue(MESH_parser,'elements',0);
      addParamValue(MESH_parser,'numNodes',0);
else
      addParameter(MESH_parser,'nodes',0);
      addParameter(MESH_parser,'boundaries',0);
      addParameter(MESH_parser,'elements',0);
      addParameter(MESH_parser,'numNodes',0);
end

parse(MESH_parser,MESH);

if ~isempty(MESH_parser.UsingDefaults)
    msg = sprintf('MESH data structure does not contain the following fields: \n');
    for i = 1 : length(MESH_parser.UsingDefaults)
        msg = [msg, ' ', MESH_parser.UsingDefaults{i}];
    end
    error(msg)
end


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
        Dirichlet_side          = [Dirichlet_side,find(MESH.boundaries(5,:) == type_Dirichlet(k))];
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
        Neumann_side = [Neumann_side,find(MESH.boundaries(5,:) == type_Neumann(k))];
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
        Robin_side = [Robin_side,find(MESH.boundaries(5,:) == type_Robin(k))];
    end
    MESH.Robin_side = unique(Robin_side);
    
else
    MESH.Robin_side = [];
end



end