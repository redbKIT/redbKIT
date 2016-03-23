function [ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA )
%BUILDMESH generates FE_SPACE struct

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% Fill MESH data structure
MESH.dim         = dim;
MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements    = elements;
MESH.numVertices = size(vertices,2);

%% Build higher order (P2 or P3) mesh if required
if ~strcmp(fem,'P1')
    fprintf('\n Generating P2 mesh ... ')
    time_mesh = tic;
    [MESH.elements, MESH.nodes, MESH.boundaries] = ...
        feval(['P1to',fem,'mesh',num2str(dim),'D'],elements, vertices, boundaries);
    time_mesh = toc(time_mesh);
    fprintf('done in %f s\n', time_mesh)
else
    MESH.nodes = vertices;
end


%% Update Mesh data with BC information and geometrical maps
[~,numBoundaryDof]           = select(fem, dim);
MESH.numNodes                = size(MESH.nodes,2);
MESH.numElem                 = size(MESH.elements,2);
MESH.numBoundaryDof          = numBoundaryDof;

% Update MESH with BC information
[MESH]         = BC_info(MESH, DATA);

% Compute geometrical map (ref to physical elements) information
[MESH.jac, MESH.invjac, MESH.h] = geotrasf(dim, MESH.vertices, MESH.elements);   
    
% Compute quadrature nodes and weights on the reference element
[quad_nodes]  = quadrature(dim, quad_order);

% Evaluate P1 geometrical mapping basis functions in the quad points
[MESH.chi]                  =  fem_basis(dim, 'P1', quad_nodes);

end

