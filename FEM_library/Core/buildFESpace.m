function [ FE_SPACE ] = buildFESpace( MESH, fem, numComponents, quad_order )
%BUILDFESPACE generates FE_SPACE struct

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

[numElemDof,numBoundaryDof]  = select(fem, MESH.dim);
[quad_nodes, quad_weights]   = quadrature(MESH.dim, quad_order);

FE_SPACE.dim              = MESH.dim;
FE_SPACE.fem              = fem;
FE_SPACE.numComponents    = numComponents;
FE_SPACE.numElemDof       = numElemDof;
FE_SPACE.numBoundaryDof   = numBoundaryDof;

if strcmp( fem , 'P1')
    FE_SPACE.numDof           = numComponents * MESH.numVertices;
    FE_SPACE.numDofScalar     = MESH.numVertices;
else
    FE_SPACE.numDof           = numComponents * MESH.numNodes;
    FE_SPACE.numDofScalar     = MESH.numNodes;
end

% Store quadrature nodes and weights on the reference element
FE_SPACE.quad_order    = quad_order;
FE_SPACE.quad_nodes    = quad_nodes;
FE_SPACE.quad_weights  = quad_weights;
FE_SPACE.numQuadNodes  = length(FE_SPACE.quad_nodes);

% Evaluate basis functions in the quad points on the reference element
[FE_SPACE.phi, FE_SPACE.dphi_ref] =  fem_basis(FE_SPACE.dim, FE_SPACE.fem, FE_SPACE.quad_nodes);

end

