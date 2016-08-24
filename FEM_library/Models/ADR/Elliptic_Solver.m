function [u, FE_SPACE, MESH, DATA, errorL2, errorH1] = Elliptic_Solver(dim, elements, vertices, boundaries, fem, data_file, param, vtk_filename)
%ELLIPTIC_SOLVER diffusion-transport-reaction finite element solver
%
%   [U, FE_SPACE, MESH, DATA, ERRORL2, ERRORH1] = ...
%    ELLIPTIC2D_SOLVER(DIM, ELEMENTS, VERTICES, BOUNDARIES, FEM, DATA_FILE, 
%                      PARAM, VTK_FILENAME)
%
%   Inputs:
%     DIM: space dimension, either 2 or 3
%     ELEMENTS, VERTICES, BOUNDARIES: mesh information
%     FEM: string 'P1' or 'P2'
%     DATA_FILE: name of the file defining the problem data and
%          boundary conditions.
%     PARAM: vector of parameters possibly used in the data_file; 
%         if not provided, the PARAM vector is set to the empty vector.
%     VTK_FILENAME: string containing the filename for exporting the
%         solution in the VTK File Format. If not provided or empty, the
%         solution is not exported to vtk.
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

if nargin < 7
    param = [];
end

if nargin < 8
    vtk_filename = [];
end

%% Read problem parameters and BCs from data_file
DATA   = read_DataFile(data_file, dim, param);
DATA.param = param;

use_SUPG = false;
if isfield(DATA, 'Stabilization')
    if strcmp( DATA.Stabilization, 'SUPG' )
        if strcmp(fem, 'P1')
            use_SUPG = true;
        else
            warning('SUPG Stabilization available only for P1 FEM');
        end
    end
end

%% Set quad_order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

%% Create and fill the MESH data structure
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA );

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, 1, quad_order );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf('-------------------------------------------\n');

%% Generate Domain Decomposition (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

if isfield(DATA.Preconditioner, 'type') && strcmp( DATA.Preconditioner.type, 'AdditiveSchwarz')
    
    if isfield(DATA.Preconditioner, 'coarse_level')
        if ~strcmp( DATA.Preconditioner.coarse_level, 'None')
            R = ADR_overlapping_DD(MESH, DATA.Preconditioner.num_subdomains, ...
                DATA.Preconditioner.overlap_level, DATA.Preconditioner.coarse_num_aggregates);
        else
            R = ADR_overlapping_DD(MESH, DATA.Preconditioner.num_subdomains, DATA.Preconditioner.overlap_level);
        end
    else
        R = ADR_overlapping_DD(MESH, DATA.Preconditioner.num_subdomains, DATA.Preconditioner.overlap_level);
    end
    
    Precon.SetRestrictions( R );
    clear R;
end

%% Assemble matrix and right-hand side
fprintf('\n Assembling ... ');
t_assembly = tic;
[A, F]  =  ADR_Assembler(MESH, DATA, FE_SPACE);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

if use_SUPG
    fprintf('\n Assembling SUPG Terms ... ');
    t_assembly = tic;
    [A_SUPG, F_SUPG] = ADR_Assembler(MESH, DATA, FE_SPACE, [], [], [], [], [], 'SUPG');
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    A = A + A_SUPG;
    F = F + F_SUPG;
end

%% Apply boundary conditions
fprintf('\n Apply boundary conditions ');
[A_in, F_in, u_D]   =  ADR_ApplyBC(A, F, FE_SPACE, MESH, DATA);

%% Solve
LinSolver = LinearSolver( DATA.LinearSolver );
u                         = zeros(MESH.numNodes,1);

fprintf('\n Solve Au = f ... ');
Precon.Build( A_in );
fprintf('\n       **  time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
LinSolver.SetPreconditioner( Precon );
u(MESH.internal_dof)      = LinSolver.Solve( A_in, F_in );
fprintf('\n       ** time to solve the linear system in %3.3f s \n\n', LinSolver.GetSolveTime());

u(MESH.Dirichlet_dof)     = u_D;

%% Export to VTK
if ~isempty(vtk_filename)
    ADR_export_solution(MESH.dim, u(1:MESH.numVertices), MESH.vertices, MESH.elements, vtk_filename);
end

%% Compute L2 and H1 errors
errorL2 = [];
errorH1 = [];

if nargout == 5
    [errorL2] = FEM_error(u, MESH, DATA, FE_SPACE);
    fprintf(' L2-error : %1.3e\n', errorL2);
elseif nargout == 6
    [errorL2,errorH1] = FEM_error(u, MESH, DATA, FE_SPACE);
    fprintf(' L2-error : %1.3e H1-error : %1.3e\n',errorL2, errorH1);
end

%% Store matrix and rhs into FE_SPACE struct
FE_SPACE.A_in = A_in;
FE_SPACE.F_in = F_in;

return
