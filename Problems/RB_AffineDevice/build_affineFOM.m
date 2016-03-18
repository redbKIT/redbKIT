function [ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, data_file )
%   For reference, see Section 8.3 and 9.1 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% For the moment HARDCODED inputs
FOM.Qa = 6;
FOM.Qf = 1;
FOM.P  = 3;
FOM.mu_min = [ 0.0      1      2     ];
FOM.mu_max = [ 0.6      20     50    ];

%% Read problem parameters and BCs from data_file
DATA   = read_DataFile(data_file);
if nargin < 6
    DATA.param = [];
else
    DATA.param = param;
end
t      = [];


%% Fill MESH data structure
dim              = 2;
MESH.dim         = dim;
MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements    = elements;
MESH.numVertices = size(vertices,2);

%% Build higher order (P2 or P3) mesh if required
if ~strcmp(fem,'P1')
    [MESH.elements, MESH.nodes, MESH.boundaries] = ...
        feval(strcat('P1to',fem,'mesh','2D'),elements,vertices, boundaries);
else
    MESH.nodes = vertices;
end


%% Update Mesh data with BC information and geometrical maps
[numElemDof,numBoundaryDof]  = select(fem, dim);
MESH.numNodes                = size(MESH.nodes,2);
MESH.numElem                 = size(MESH.elements,2);
MESH.numBoundaryDof          = numBoundaryDof;

% Update MESH with BC information
[MESH]         = BC_info(MESH, DATA);

% Compute geometrical map (ref to physical elements) information
[MESH.jac, MESH.invjac, MESH.h] = geotrasf(MESH.dim, MESH.vertices, MESH.elements);   

% Compute quadrature nodes and weights on the reference element
quad_order                  = 4;  
[quad_nodes, quad_weights]  = quadrature(dim, quad_order);

% Evaluate P1 geometrical mapping basis functions in the quad points
[MESH.chi]                  =  fem_basis(dim, 'P1', quad_nodes);

%% Create and fill the FE_SPACE data structure
FE_SPACE.fem              = fem;
FE_SPACE.numDof           = length(MESH.internal_dof);
FE_SPACE.numElemDof       = numElemDof;
FE_SPACE.numBoundaryDof   = numBoundaryDof;

% Store quadrature nodes and weights on the reference element
FE_SPACE.quad_order    = quad_order;
FE_SPACE.quad_nodes    = quad_nodes;
FE_SPACE.quad_weights  = quad_weights;
FE_SPACE.numQuadNodes  = length(FE_SPACE.quad_nodes);

% Evaluate basis functions in the quad points on the reference element
[FE_SPACE.phi, FE_SPACE.dphi_ref]  =  ...
    fem_basis(dim, FE_SPACE.fem, FE_SPACE.quad_nodes);


fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',FE_SPACE.numDof);
fprintf('-------------------------------------------\n');

%% Fill FOM structures
FOM.MESH     = MESH;
FOM.FE_SPACE = FE_SPACE;
FOM.DATA     = DATA;

%% Assemble affine matrices and rhs vectors

% initialize
for i = 1 : FOM.Qa
    FOM.Aq{i} = sparse(FOM.FE_SPACE.numDof, FOM.FE_SPACE.numDof);
end

for i = 1 : FOM.Qf
    FOM.Fq{i} = sparse(FOM.FE_SPACE.numDof, 1);
end

% A_1: diffusion
A_1       =  Assembler_2D(MESH, DATA, FE_SPACE, 'diffusion', [1 1], [], 1);
FOM.Aq{1} =  A_1(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_2: diffusion
A_2       =  Assembler_2D(MESH, DATA, FE_SPACE, 'diffusion', [2 2], [], 1);
FOM.Aq{2} =  A_2(FOM.MESH.internal_dof,FOM.MESH.internal_dof);


% A_3: transport
A_3       =  Assembler_2D(MESH, DATA, FE_SPACE, 'transport', [], 2, 1);
FOM.Aq{3} =  A_3(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_4: diffusion
A_4       =  Assembler_2D(MESH, DATA, FE_SPACE, 'diffusion', [], [], 2);
FOM.Aq{4} =  A_4(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_5: diffusion
A_5       =  Assembler_2D(MESH, DATA, FE_SPACE, 'diffusion', [], [], 3);
FOM.Aq{5} =  A_5(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_6: diffusion
A_6       =  Assembler_2D(MESH, DATA, FE_SPACE, 'diffusion', [], [], 4);
FOM.Aq{6} =  A_6(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% F_1: diffusion
[~, F_1]       =  Assembler_2D(MESH, DATA, FE_SPACE, 'source', [], [], 2);
[~, F_1, u_D]  =  ApplyBC([], F_1, FE_SPACE, MESH, DATA);
FOM.Fq{1}      =  F_1;

FOM.u_D        =  @(x,mu)(FOM.DATA.bcDir(x(1,:),x(2,:),[],mu));


%% Compute Xnorm
DATA.diffusion =  @(x,y,t,param)(1+0.*x.*y);
X              =  Assembler_2D(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
DATA.reaction  =  @(x,y,t,param)(1+0.*x.*y);
M              =  Assembler_2D(MESH, DATA, FE_SPACE, 'reaction');

FOM.Xnorm =  X(FOM.MESH.internal_dof,FOM.MESH.internal_dof);
end
