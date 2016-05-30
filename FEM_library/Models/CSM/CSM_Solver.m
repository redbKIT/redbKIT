function [u, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, data_file, param, vtk_filename)
%CSM_SOLVER Static Structural Finite Element Solver

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 6
    error('Missing input arguments. Please type help CSM_Solver')
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
DATA   = CSM_read_DataFile(data_file, dim, param);
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

%% Generate Domain Decomposition (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

if isfield(DATA.Preconditioner, 'type') && strcmp( DATA.Preconditioner.type, 'AdditiveSchwarz')
    R      = CSM_overlapping_DD(MESH, DATA.Preconditioner.num_subdomains,  DATA.Preconditioner.overlap_level);
    Precon.SetRestrictions( R );
end

%% Newton Method

tol        = DATA.NonLinearSolver.tol;
resRelNorm = tol + 1;
incrNorm   = tol + 1;
maxIter    = DATA.NonLinearSolver.maxit;
k          = 1;

[~, ~, u_D]   =  CSM_ApplyBC([], [], FE_SPACE, MESH, DATA);
dU             = zeros(MESH.numNodes*MESH.dim,1);
U_k            = zeros(MESH.numNodes*MESH.dim,1);
U_k(MESH.Dirichlet_dof) = u_D;

% Assemble matrix and right-hand side
fprintf('\n -- Assembling external Forces... ');
t_assembly = tic;
F_ext      = CSM_Assembler('external_forces', MESH, DATA, FE_SPACE);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

fprintf('\n -- Assembling internal Forces... ');
t_assembly = tic;
[F_in, dF_in]  =  CSM_Assembler('internal_forces', MESH, DATA, FE_SPACE, U_k);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

Residual = F_in - F_ext;
% Apply boundary conditions
fprintf('\n -- Apply boundary conditions ... ');
t_assembly = tic;
[A, b]   =  CSM_ApplyBC(dF_in, -Residual, FE_SPACE, MESH, DATA, [], 1);
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
    
    % update solution
    U_k        = U_k + dU;
    incrNorm   = norm(dU)/norm(U_k);
    
    % Assemble matrix and right-hand side
    fprintf('\n   -- Assembling internal forces... ');
    t_assembly = tic;
    [F_in, dF_in]  =  CSM_Assembler('internal_forces', MESH, DATA, FE_SPACE, full(U_k));
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    Residual   = F_in - F_ext;
    
    % Apply boundary conditions
    fprintf('\n   -- Apply boundary conditions ... ');
    t_assembly = tic;
    [A, b]   =  CSM_ApplyBC(dF_in, -Residual, FE_SPACE, MESH, DATA, [], 1);
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    resRelNorm = norm(b) / res0Norm;
    
    fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k,full(incrNorm), full(norm(resRelNorm)));
    k = k + 1;
    
end

u = U_k;

%% Export Displacement to VTK
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename);
end

%% Compute Von Mises Stress
if DATA.Output.ComputeVonMisesStress
    fprintf('\n   -- Compute Element Stresses... ');
    t_assembly = tic;
    [Sigma]  =  CSM_Assembler('stress', MESH, DATA, FE_SPACE, full(U_k));
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    Sigma_VM = sqrt( (Sigma(:,1) - Sigma(:,5)).^2 + (Sigma(:,5) - Sigma(:,9)).^2 + (Sigma(:,9) - Sigma(:,1)).^2 );
    CSM_export_VonMisesStress(MESH.dim, Sigma_VM, MESH.vertices, MESH.elements, [vtk_filename, '_VMstress']);
end
 
return
