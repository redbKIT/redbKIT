function [ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA, ...
    model, rings, reduced_elements, reduced_boundaries )
%BUILDMESH generates MESH struct

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 8
    model = [];
end
%% Fill MESH data structure
MESH.dim         = dim;
MESH.fem         = fem;
MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements    = elements;
MESH.numVertices = size(vertices,2);
MESH.quad_order  = quad_order;

if nargin > 8
    MESH.rings    = rings;
end

%% Build higher order (P2 or P3) mesh if required
if ~strcmp(fem,'P1')
    fprintf('\n Generating %s mesh ... ', fem)
    time_mesh = tic;
    if nargin > 8
        [MESH.elements, MESH.nodes, MESH.boundaries, MESH.rings] = ...
        feval(['P1to',fem,'mesh',num2str(dim),'D'],elements, vertices, boundaries, rings);
    else
        [MESH.elements, MESH.nodes, MESH.boundaries] = ...
        feval(['P1to',fem,'mesh',num2str(dim),'D'],elements, vertices, boundaries);
    end
    time_mesh = toc(time_mesh);
    fprintf('done in %f s\n', time_mesh)
else
    MESH.nodes = vertices;
end

if nargin > 9
    MESH.elements = MESH.elements(:, reduced_elements);
end
    
% I cannot just restict boundaries to the reduced ones, otherwise essential
% Dofs and internal numbering are not correct.
%
% if nargin > 10    
%     MESH.boundaries = MESH.boundaries(:, reduced_boundaries);
% end

%% Update Mesh data with geometrical maps
[MESH.numElemDof,MESH.numBoundaryDof,MESH.numRingsDof]    = select(fem, dim);
MESH.numNodes                = size(MESH.nodes,2);
MESH.numElem                 = size(MESH.elements,2);

% Compute geometrical map (ref to physical elements) information
[MESH.jac, MESH.invjac, MESH.h] = geotrasf(dim, MESH.vertices, MESH.elements);   
    
% Compute quadrature nodes and weights on the reference element
[quad_nodes]  = quadrature(dim, quad_order);

% Evaluate P1 geometrical mapping basis functions in the quad points
[MESH.chi]                  =  fem_basis(dim, 'P1', quad_nodes);


%% Generate mesh normals
if strcmp( model, 'CSM') || strcmp( model, 'CFD')
    fprintf('\n Generating mesh normals ... ')
    time_mesh = tic;
    switch dim
        case 2
            [MESH.Normal_Faces] = ComputeSurfaceNormals2D(MESH.boundaries(1:2,:),MESH.vertices(1:2,:),elements(1:3,:));
            
        case 3
            [MESH.Normal_Faces] = ...
                ComputeSurfaceNormals3D(MESH.boundaries(1:3,:),MESH.vertices(1:3,:), elements(1:4,:));
    end
    time_mesh = toc(time_mesh);
    fprintf('done in %f s\n', time_mesh)
end

%% Update Mesh data with BC information
if nargin >= 7 && ~isempty(DATA)
    % Update MESH with BC information
    [MESH]         = BC_info(MESH, DATA, model);
end

end

