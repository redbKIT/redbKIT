function [u, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, data_file, param)
%CSM_SOLVER Static Structural Finite Element Solver
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
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 6
    error('Missing input arguments. Please type help CSM_Solver')
end

if isempty(data_file)
    error('Missing data_file')
end

%% Read problem parameters and BCs from data_file
DATA   = CSM_read_DataFile(data_file);
if nargin < 7
    DATA.param = [];
else
    DATA.param = param;
end
t      = [];

%% Set quad_order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

%% Create and fill the MESH data structure
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA, 'CSM' );

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, dim, quad_order );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',length(MESH.internal_dof));
fprintf('-------------------------------------------\n');


%% Assemble matrix and right-hand side
fprintf('\n >> Assembling ... ');
t_assembly = tic;
[F, A]  =  CSM_Assembler('all', MESH, DATA, FE_SPACE);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);


%% Apply boundary conditions
fprintf('\n >> Apply boundary conditions ... ');
t_assembly = tic;
[A_in, F_in, u_D]   =  CSM_ApplyBC(A, -F, FE_SPACE, MESH, DATA);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

%% Solve
fprintf('\n >> Solve Au = f ... ');
t_solve = tic;
u                         = zeros(MESH.numNodes*MESH.dim,1);
u(MESH.internal_dof)      = A_in \ F_in;
u(MESH.Dirichlet_dof)     = u_D;
t_solve = toc(t_solve);
fprintf('done in %3.3f s \n', t_solve);

STR_export_solution(MESH.dim, u, MESH.vertices, MESH.elements, MESH.numVertices, 'SOL_ELA');

%% Store matrix and rhs into FE_SPACE struct
%FE_SPACE.A_in = A_in;
%FE_SPACE.F_in = F_in;

return
