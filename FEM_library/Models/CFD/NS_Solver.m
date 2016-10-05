function [u, MESH, DATA] = NS_Solver(dim, elements, vertices, boundaries, fem, data_file, param, vtk_filename)
%NS_SOLVER steady Navier-Stokes Equations solver

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 6
    error('Missing input arguments. Please type help NSsteadySolver')
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
DATA   = CFD_read_DataFile(data_file, dim, param);
if nargin < 7
    DATA.param = [];
else
    DATA.param = param;
end
t      = [];

use_SUPG = false;
if isfield(DATA, 'Stabilization')
    if strcmp( DATA.Stabilization, 'SUPG' ) && strcmp(fem{1}, 'P1')
        use_SUPG = true;
    end
end

%% Set quad_order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

%% Create and fill the MESH data structure
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem{1}, quad_order, DATA, 'CFD' );

%% Create and fill the FE_SPACE data structure
[ FE_SPACE_v ] = buildFESpace( MESH, fem{1}, dim, quad_order );
[ FE_SPACE_p ] = buildFESpace( MESH, fem{2}, 1, quad_order );

MESH.internal_dof_c{MESH.dim+1} = 1:FE_SPACE_p.numDof;
        
totSize = FE_SPACE_v.numDof + FE_SPACE_p.numDof;

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',length(MESH.internal_dof));
fprintf('-------------------------------------------\n');

%% Generate Domain Decomposition (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

if isfield(DATA.Preconditioner, 'type') && strcmp( DATA.Preconditioner.type, 'AdditiveSchwarz')
    R      = CFD_overlapping_DD(MESH, FE_SPACE_v, FE_SPACE_p, DATA.Preconditioner.num_subdomains,  DATA.Preconditioner.overlap_level);
    Precon.SetRestrictions( R );
end

%% Create Fluid Assembler Object
FluidModel = CFD_Assembler( MESH, DATA, FE_SPACE_v, FE_SPACE_p );

%% Assemble Constant Terms
fprintf('\n   -- Assembling Stokes terms... ');
t_assembly = tic;
[A_Stokes] = FluidModel.compute_Stokes_matrix();
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);


%% Nonlinear Iterations
tol        = DATA.NonLinearSolver.tol;
resRelNorm = tol + 1;
incrNorm   = tol + 1;
maxIter    = DATA.NonLinearSolver.maxit;
k          = 1;

[~, ~, u_D]   =  CFD_ApplyBC([], [], FE_SPACE_v, FE_SPACE_p, MESH, DATA);
dU             = zeros(totSize,1);
U_k            = zeros(totSize,1);
U_k(MESH.Dirichlet_dof) = u_D;

% Assemble matrix and right-hand side
fprintf('\n   -- Assembling Convective terms... ');
t_assembly = tic;
[C1, C2] = FluidModel.compute_convective_matrix( U_k );
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

Residual = A_Stokes * U_k + C1 * U_k;
Jacobian = A_Stokes + C1 + C2;

if use_SUPG
    fprintf('\n   -- Assembling SUPG Terms ... ');
    t_assembly = tic;
    [A_SUPG, F_SUPG] = FluidModel.compute_SUPG_implicitSteady( U_k );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    Jacobian    = Jacobian + A_SUPG;
    Residual    = Residual + F_SUPG;
end

% Apply boundary conditions
fprintf('\n   -- Apply boundary conditions ... ');
t_assembly = tic;
[A, b]   =  CFD_ApplyBC(Jacobian, -Residual, FE_SPACE_v, FE_SPACE_p, MESH, DATA, [], 1);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

res0Norm = norm(b);

LinSolver = LinearSolver( DATA.LinearSolver );

fprintf('\n============ Start Newton Iterations ============\n\n');
while (k <= maxIter && incrNorm > tol && resRelNorm > tol)
            
    % Solve
    fprintf('\n   -- Solve J x = -R ... ');    
    Precon.Build( A );
    fprintf('\n        time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
    LinSolver.SetPreconditioner( Precon );
    dU(MESH.internal_dof) = LinSolver.Solve( A, b );
    fprintf('\n        time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());

    U_k        = U_k + dU;
    incrNorm   = norm(dU)/norm(U_k);
    
    % Assemble matrix and right-hand side
    fprintf('\n   -- Assembling Convective terms... ');
    t_assembly = tic;
    [C1, C2] = FluidModel.compute_convective_matrix( U_k );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    Residual = A_Stokes * U_k + C1 * U_k;
    Jacobian = A_Stokes + C1 + C2;
    
    if use_SUPG
        fprintf('\n   -- Assembling SUPG Terms ... ');
        t_assembly = tic;
        [A_SUPG, F_SUPG] = FluidModel.compute_SUPG_implicitSteady( U_k );
        t_assembly = toc(t_assembly);
        fprintf('done in %3.3f s\n', t_assembly);
        
        Jacobian    = Jacobian + A_SUPG;
        Residual    = Residual + F_SUPG;
    end
    
    % Apply boundary conditions
    fprintf('\n   -- Apply boundary conditions ... ');
    t_assembly = tic;
    [A, b]   =  CFD_ApplyBC(Jacobian, -Residual, FE_SPACE_v, FE_SPACE_p, MESH, DATA, [], 1);
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    resRelNorm = norm(b) / res0Norm;
    
    fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k,full(incrNorm), full(norm(resRelNorm)));
    k = k + 1;
    
end

u = U_k;

%% Export to VTK
if ~isempty(vtk_filename)
    CFD_export_solution(MESH.dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename);
end
 
return
