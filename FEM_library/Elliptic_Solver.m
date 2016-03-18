function [u, FE_SPACE, MESH, DATA, errorL2, errorH1] = Elliptic_Solver(dim, elements, vertices, boundaries, fem, data_file, param)
%ELLIPTIC2D_SOLVER 2D diffusion-transport-reaction finite element solver
%
%   [U, FE_SPACE, MESH, DATA, ERRORL2, ERRORH1] = ...
%    ELLIPTIC2D_SOLVER(ELEMENTS, VERTICES, BOUNDARIES, FEM, DATA_FILE, PARAM)
%
%   Inputs:
%     ELEMENTS, VERTICES, BOUNDARIES: mesh information
%     FEM: string 'P1' or 'P2'
%     DATA_FILE: name of the file defining the problem data and
%          boundary conditions.
%     PARAM: vector of parameters possibly used in the data_file; 
%         if not provided, the PARAM vector is set to the empty vector.
%
%   Outputs:
%     U: problem solution
%     ERRORL2: L2-error between the numerical solution and the exact one 
%        (provided by the user in the data_file)
%     ERRORH1: H1-error between the numerical solution and the exact one 
%        (provided by the user in the data_file)
%     FE_SPACE: struct containing Finite Element Space information
%     MESH: struct containing mesh information
%     DATA: struct containing problem data

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 6
    error('Missing input arguments. Please type help Elliptic_Solver')
end

if isempty(data_file)
    error('Missing data_file')
end

%% Read problem parameters and BCs from data_file
DATA   = read_DataFile(data_file);
if nargin < 7
    DATA.param = [];
else
    DATA.param = param;
end
t      = [];

%% Fill MESH data structure
MESH.dim         = dim;
MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements    = elements;
MESH.numVertices = size(vertices,2);

%% Build higher order (P2 or P3) mesh if required
if ~strcmp(fem,'P1')
    fprintf('\n Generate P2 mesh\n')
    [MESH.elements, MESH.nodes, MESH.boundaries] = ...
        feval(['P1to',fem,'mesh',num2str(dim),'D'],elements, vertices, boundaries);
    fprintf('\n Done\n')
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
[MESH.jac, MESH.invjac, MESH.h] = geotrasf(dim, MESH.vertices, MESH.elements);   

% Compute quadrature nodes and weights on the reference element
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end
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
% [FE_SPACE.phi, FE_SPACE.dcsiphi, FE_SPACE.detaphi]  =  ...
%     fem_basis(FE_SPACE.fem, FE_SPACE.quad_nodes);

[FE_SPACE.phi, FE_SPACE.dphi_ref]  =  ...
    fem_basis(dim, FE_SPACE.fem, FE_SPACE.quad_nodes);


fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf('-------------------------------------------\n');


%% Assemble matrix and right-hand side
fprintf('\n Assembling ... ');
t_assembly = tic;
switch dim
    case 2
        [A, F]  =  Assembler_2D(MESH, DATA, FE_SPACE);
    case 3
        [A, F]  =  Assembler_3D(MESH, DATA, FE_SPACE);     
end
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);


%% Apply boundary conditions
fprintf('\n Apply boundary conditions ');
[A_in, F_in, u_D]   =  ApplyBC(A, F, FE_SPACE, MESH, DATA);


%% Solve
fprintf('\n Solve Au = f ... ');
t_solve = tic;
u                         = zeros(MESH.numNodes,1);
u(MESH.internal_dof)      = A_in \ F_in;
u(MESH.Dirichlet_dof)     = u_D;t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);


%% Compute L2 and H1 errors
errorL2 = [];
errorH1 = [];

if nargout == 5
    [errorL2] = FEM_error(u, MESH, DATA, FE_SPACE);
    fprintf(' L2-error : %e\n', errorL2);
elseif nargout == 6
    [errorL2,errorH1] = FEM_error(u, MESH, DATA, FE_SPACE);
    fprintf(' L2-error : %e H1-error : %e\n',errorL2, errorH1);
end

%% Store matrix and rhs into FE_SPACE struct
FE_SPACE.A_in = A_in;
FE_SPACE.F_in = F_in;

return
