function [u, FE_SPACE, MESH, DATA] = Spinning_Solver(dim, elements, vertices, boundaries, fem, data_file, param, vtk_filename)
%CSMT_SOLVER Dynamic Structural Finite Element Solver

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

%% Gather Time Setting
t0        = DATA.time.t0;
dt        = DATA.time.dt;
tf        = DATA.time.tf;
t         = DATA.time.t0;
k_t       = 0;

%% Create and fill the MESH data structure
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA, 'CSM' );


u0  = [];
du0 = [];
for k = 1 : MESH.dim
    switch dim
        case 2
            u0  = [u0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), t0, param )'];
            du0 = [du0; DATA.du0{k}( MESH.nodes(1,:), MESH.nodes(2,:), t0, param )'];
            
        case 3
            u0  = [u0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];
            du0 = [du0; DATA.du0{k}( MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];
    end
end
d2u0 = zeros(MESH.numNodes*MESH.dim,1);
u = u0;

indx = find(isnan(du0));
du0(indx) = 0;

alpha = pi/3;
Rot_X_Matrix = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
dU_n_tmp     = [du0(1:MESH.numNodes)'; du0(1+MESH.numNodes:2*MESH.numNodes)'; du0(2*MESH.numNodes+1:end)'];
dU_n_tmp     = Rot_X_Matrix  * dU_n_tmp;
du0          = [dU_n_tmp(1,:)'; dU_n_tmp(2,:)'; dU_n_tmp(3,:)'];
d2u0 = zeros(MESH.numNodes*MESH.dim,1);

%def_vert = Rot_X_Matrix * vertices;
%[ MESH ] = buildMESH( dim, elements, def_vert, boundaries, fem, quad_order, DATA, 'CSM' );

%================= Only for Top Spinning ==================================
indexApex = find(ismember(MESH.nodes(1:3,:)',[0.0 0.0 0.0],'rows'));
MESH.Dirichlet_dof_c{1} = indexApex;
MESH.Dirichlet_dof_c{2} = indexApex;
MESH.Dirichlet_dof_c{3} = indexApex;

MESH.Dirichlet_dof = [MESH.Dirichlet_dof_c{1} MESH.numNodes+MESH.Dirichlet_dof_c{2} 2*MESH.numNodes+MESH.Dirichlet_dof_c{3}];

MESH.internal_dof = [];
for k = 1 : MESH.dim
    MESH.internal_dof_c{k}   = setdiff([1:MESH.numNodes],MESH.Dirichlet_dof_c{k});
    MESH.internal_dof        = [MESH.internal_dof (k-1)*MESH.numNodes+MESH.internal_dof_c{k}];
end
%==========================================================================

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, dim, quad_order );

TimeAdvance = Newmark_TimeAdvance( DATA.time.beta, DATA.time.gamma, dt );

CSM_export_solution(MESH.dim, u0, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, 0);
CSM_export_solution(MESH.dim, du0, MESH.vertices, MESH.elements, MESH.numNodes, 'Figures/Initial_velocity');

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',length(MESH.internal_dof));
fprintf(' * Number of timesteps =  %d\n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

%% Generate Domain Decomposition (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

if isfield(DATA.Preconditioner, 'type') && strcmp( DATA.Preconditioner.type, 'AdditiveSchwarz')
    R      = CSM_overlapping_DD(MESH, DATA.Preconditioner.num_subdomains,  DATA.Preconditioner.overlap_level);
    Precon.SetRestrictions( R );
end

%% Assemble mass matrix
fprintf('\n Assembling mass matrix... ');
t_assembly = tic;
M    =  CSM_Assembler('mass', MESH, DATA, FE_SPACE);
M    =  M * DATA.Density;
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

LinSolver = LinearSolver( DATA.LinearSolver );

%d2u0 = (1 / DATA.Density * M) \ du0;

fprintf('\n -- Assembling external Forces... ');
t_assembly = tic;
F_ext      = CSM_Assembler('external_forces', MESH, DATA, FE_SPACE, [], t0);%M*(-9.81*ones(size(M,1),1));%
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

fprintf('\n -- Assembling internal Forces... ');
t_assembly = tic;
[F_in]  =  CSM_Assembler('internal_forces', MESH, DATA, FE_SPACE, u0);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

%d2u0 = (1 / DATA.Density * M) \ (F_ext - F_in);

TimeAdvance.Initialize( u0, du0, d2u0 );
Coef_Mass = TimeAdvance.MassCoefficient( );


%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );

    % Newton Method
    tol        = DATA.NonLinearSolver.tol;
    resRelNorm = tol + 1;
    incrNorm   = tol + 1;
    maxIter    = DATA.NonLinearSolver.maxit;
    k          = 1;
    
    [~, ~, u_D]   =  CSM_ApplyBC([], [], FE_SPACE, MESH, DATA, t);
    dU             = zeros(MESH.numNodes*MESH.dim,1);
    U_k            = u;
    U_k(MESH.Dirichlet_dof) = u_D;
    
    Csi = TimeAdvance.RhsContribute( );
    
    % Assemble matrix and right-hand side
    fprintf('\n -- Assembling external Forces... ');
    t_assembly = tic;
    F_ext      = CSM_Assembler('external_forces', MESH, DATA, FE_SPACE, [], t);%M*(-9.81*ones(size(M,1),1));%
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    F_ext = F_ext + M * Csi;
    
    fprintf('\n -- Assembling internal Forces... ');
    t_assembly = tic;
    [F_in, dF_in]  =  CSM_Assembler('internal_forces', MESH, DATA, FE_SPACE, U_k);
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    Residual  = Coef_Mass * M * U_k +  F_in - F_ext;
    Jacobian  = Coef_Mass * M + dF_in;
    
    % Apply boundary conditions
    fprintf('\n -- Apply boundary conditions ... ');
    t_assembly = tic;
    [A, b]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE, MESH, DATA, t, 1);
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
        
    res0Norm = norm(b);
        
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
        
        Residual  = Coef_Mass * M * U_k +  F_in - F_ext;
        Jacobian  = Coef_Mass * M + dF_in;
        
        % Apply boundary conditions
        fprintf('\n   -- Apply boundary conditions ... ');
        t_assembly = tic;
        [A, b]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE, MESH, DATA, t, 1);
        t_assembly = toc(t_assembly);
        fprintf('done in %3.3f s\n', t_assembly);
        
        resRelNorm = norm(b) / res0Norm;
        
        fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k,full(incrNorm), full(norm(resRelNorm)));
        k = k + 1;
        
    end
    
    u = U_k;
    
    %% Export to VTK
    if ~isempty(vtk_filename) && mod(k_t,2)
        CSM_export_solution(MESH.dim, u, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, k_t);
    end
    
    TimeAdvance.Update( u );
    
    iter_time = toc(iter_time);
    fprintf('\n-------------- Iteration time: %3.2f s -----------------',iter_time);
    
end

fprintf('\n************************************************************************* \n');

return
