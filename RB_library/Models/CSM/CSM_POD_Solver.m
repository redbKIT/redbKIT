function [u, FE_SPACE, MESH, DATA] = CSM_POD_Solver(dim, elements, vertices, boundaries, fem, data_file, param, vtk_filename, InitialGuess, Training_Options, V_POD)
%CSM_POD_SOLVER Static Structural POD-Galerkin Solver
%
%   Training_Options should be a struct with fields:
%      - Training_Options.h5_filename
%      - Training_Options.InternalForces.h5_section
%      - Training_Options.ExternalForces.h5_section
%
%   V_POD is a POD basis ( i.e. matrix of size #InternalDoFs x #PODmodes )

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

if nargin < 8 || isempty( vtk_filename )
    vtk_filename = [];
end

if nargin < 9 || isempty( InitialGuess )
    InitialGuess = [];
end

if nargin < 10 || isempty( Training_Options )
    export_h5 = false;
else
    export_h5    = true;
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
fprintf(' * Number of Dofs      = %d \n',size(V_POD, 2));
fprintf('-------------------------------------------\n');

if export_h5
    Fint_h5  = HDF5_DenseMultiCVector(Training_Options.h5_filename, ...
        Training_Options.InternalForces.h5_section, length(MESH.internal_dof));
    Fext_h5  = HDF5_DenseMultiCVector(Training_Options.h5_filename, ...
        Training_Options.ExternalForces.h5_section, length(MESH.internal_dof));
end

%% Generate Domain Decomposition (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

%% Newton Method
tolNewton  = DATA.NonLinearSolver.tol;
resRelNorm = tolNewton + 1;
incrNorm   = tolNewton + 1;
maxIter    = DATA.NonLinearSolver.maxit;
k          = 1;
backtrackMaxIter   = DATA.NonLinearSolver.backtrackIter;
backtrackFactor    = DATA.NonLinearSolver.backtrackFactor;

[~, ~, u_D]   =  CSM_ApplyBC([], [], FE_SPACE, MESH, DATA);
dU             = zeros(FE_SPACE.numDof,1);

% Initial Guess
if isempty( InitialGuess )
    U_k            = zeros(FE_SPACE.numDof,1);
    U_k(MESH.Dirichlet_dof) = u_D;
else
    if size(InitialGuess,1 ) == FE_SPACE.numDof && size(InitialGuess,2 ) == 1
        U_k = InitialGuess;
    else
        error('InitialGuess should be of size FE_SPACE.numDofs x 1')
    end
end

SolidModel = CSM_Assembler( MESH, DATA, FE_SPACE );

% Assemble matrix and right-hand side
fprintf('\n -- Assembling external Forces... ');
t_assembly = tic;
F_ext      = SolidModel.compute_volumetric_forces();
[~, F_ext] = CSM_ApplyBC([], F_ext, FE_SPACE, MESH, DATA, [], 1);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

if export_h5
    Fext_h5.append( F_ext );
end

% Assemble Robin BC
fprintf('\n -- Assembling Robin BC... ');
t_assembly = tic;
A_robin = SolidModel.assemble_ElasticRobinBC();
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

fprintf('\n -- Assembling internal Forces... ');
t_assembly = tic;
F_in      = SolidModel.compute_internal_forces(U_k);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

F_in = F_in + A_robin * U_k;

if export_h5
    Fint_h5.append( F_in(MESH.internal_dof) );
end

Residual = F_in(MESH.internal_dof) - F_ext;
res0Norm = norm(full(Residual));
resNorm_old = norm(full(Residual));

LinSolver = LinearSolver( DATA.LinearSolver );

fprintf('\n============ Start Newton Iterations ============\n\n');
while (k <= maxIter && incrNorm > tolNewton && resRelNorm > tolNewton)
    
    fprintf('\n -- Assembling jacobian matrix... ');
    t_assembly = tic;
    dF_in     = SolidModel.compute_jacobian(U_k);
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    % Apply boundary conditions
    fprintf('\n -- Apply boundary conditions ... ');
    t_assembly = tic;
    A   =  CSM_ApplyBC(dF_in + A_robin, [], FE_SPACE, MESH, DATA, [], 1);
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);

    % Solve
    fprintf('\n   -- Solve J x = -R ... ');  
    A_N = V_POD' * (A * V_POD);
    Precon.Build( A_N );
    fprintf('\n        time to build the preconditioner: %3.3f s \n', Precon.GetBuildTime());
    LinSolver.SetPreconditioner( Precon );
    dU_N                  = LinSolver.Solve( A_N, -V_POD' * Residual );
    dU(MESH.internal_dof) = V_POD * dU_N;
    fprintf('\n        time to solve the linear system: %3.3f s \n', LinSolver.GetSolveTime());

    % update Zero Normal Displacement condition
    %dU = MESH.DirichletNormal_R * dU;
    
    % update solution
    U_k_tmp     = U_k + dU;
    
    % Assemble residual
    fprintf('\n   -- Assembling internal forces... ');
    t_assembly = tic;
    F_in       = SolidModel.compute_internal_forces(U_k_tmp);
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    F_in = F_in + A_robin * U_k_tmp;
    Residual = F_in(MESH.internal_dof) - F_ext;
    resRelNorm = norm(Residual) / res0Norm;

    % backtracking if needed
    alpha          = 1;
    backtrack_iter = 0;
    
    while ((norm(full(Residual)) > 2 * resNorm_old || isnan( norm(full(Residual)) )) && backtrack_iter < backtrackMaxIter )
    
        alpha = alpha * backtrackFactor;
        backtrack_iter = backtrack_iter + 1;
        
        % update solution
        U_k_tmp     = U_k + alpha * dU;
        
        % Assemble residual
        fprintf('\n    -- Backtracing alpha = %1.2e: Assembling internal forces... ', alpha);
        t_assembly = tic;
        F_in       = SolidModel.compute_internal_forces(U_k_tmp);
        t_assembly = toc(t_assembly);
        F_in     = F_in + A_robin * U_k_tmp;
        Residual = F_in(MESH.internal_dof) - F_ext;
        resRelNorm = norm(full(Residual)) / res0Norm;
        fprintf('done in %3.3f s, Residual Rel Norm = %1.2e\n', t_assembly, full(norm(resRelNorm)));
        
    end
    
    U_k = U_k_tmp;
    incrNorm   = norm(alpha*dU)/norm(U_k);
    resNorm_old = norm(Residual);
    
    if export_h5
        Fint_h5.append( F_in(MESH.internal_dof) );
    end

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
    [Sigma]  =  SolidModel.compute_stress(U_k);%CSM_Assembler_function('stress', MESH, DATA, FE_SPACE, full(U_k));
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    if MESH.dim == 2
        Sigma_VM = sqrt(  Sigma(:,1).^2 + Sigma(:,4).^2 - Sigma(:,1) .* Sigma(:,4) + 3 * Sigma(:,2).^2 );
    elseif MESH.dim == 3
        Sigma_VM = sqrt( 0.5 * ( (Sigma(:,1) - Sigma(:,5)).^2 + (Sigma(:,5) - Sigma(:,9)).^2 + (Sigma(:,9) - Sigma(:,1)).^2 + 6 * ( Sigma(:,2).^2 + Sigma(:,6).^2 + Sigma(:,7).^2 ) ) );
    end
    CSM_export_VonMisesStress(MESH.dim, Sigma_VM, MESH.vertices, MESH.elements, [vtk_filename, '_VMstress']);
end
 
return
