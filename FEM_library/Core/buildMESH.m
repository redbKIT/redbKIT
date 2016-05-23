function [ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA, model )
%BUILDMESH generates MESH struct

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% Fill MESH data structure
MESH.dim         = dim;
MESH.fem         = fem;
MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements    = elements;
MESH.numVertices = size(vertices,2);

%% Build higher order (P2 or P3) mesh if required
if ~strcmp(fem,'P1')
    fprintf('\n Generating %s mesh ... ', fem)
    time_mesh = tic;
    [MESH.elements, MESH.nodes, MESH.boundaries] = ...
        feval(['P1to',fem,'mesh',num2str(dim),'D'],elements, vertices, boundaries);
    time_mesh = toc(time_mesh);
    fprintf('done in %f s\n', time_mesh)
else
    MESH.nodes = vertices;
end


%% Update Mesh data with BC information and geometrical maps
[MESH.numElemDof,MESH.numBoundaryDof]    = select(fem, dim);
MESH.numNodes                = size(MESH.nodes,2);
MESH.numElem                 = size(MESH.elements,2);

% Compute geometrical map (ref to physical elements) information
[MESH.jac, MESH.invjac, MESH.h] = geotrasf(dim, MESH.vertices, MESH.elements);   
    
% Compute quadrature nodes and weights on the reference element
[quad_nodes]  = quadrature(dim, quad_order);

% Evaluate P1 geometrical mapping basis functions in the quad points
[MESH.chi]                  =  fem_basis(dim, 'P1', quad_nodes);

if nargin >= 7 && ~isempty(DATA)
    if nargin < 8
        model = [];
    end
    % Update MESH with BC information
    [MESH]         = BC_info(MESH, DATA, model);
end

if strcmp( model, 'CSM') || strcmp( model, 'CFD')
    fprintf('\n Generating mesh normals ... ')
    time_mesh = tic;
    switch dim
        case 2
            [MESH.Normal_Faces] = ComputeSurfaceNormals2D(MESH.boundaries(1:2,:),MESH.vertices(1:2,:),MESH.elements(1:3,:));
            
        case 3
            [MESH.Normal_Faces] = ...
                ComputeSurfaceNormals3D(MESH.boundaries(1:3,:),MESH.vertices(1:3,:), MESH.elements(1:4,:));
    end
    time_mesh = toc(time_mesh);
    fprintf('done in %f s\n', time_mesh)
end
    

end

