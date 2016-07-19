function [u, FE_SPACE, MESH, DATA] = CSMt_PODDEIM_Solver(dim, elements, vertices, boundaries, fem, data_file, ...
    param, vtk_filename, sol_history, ROM)
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

if nargin < 9 || isempty(sol_history)
    sol_history = false;
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

%% Gather Time Setting
t0        = DATA.time.t0;
dt        = DATA.time.dt;
tf        = DATA.time.tf;
t         = DATA.time.t0;
k_t       = 0;


TimeAdvance = GeneralizedAlpha_TimeAdvance( DATA.time.beta, DATA.time.gamma, DATA.time.alpha_m, DATA.time.alpha_f, dt );

u0  = [];
du0 = [];
for k = 1 : FE_SPACE.numComponents
    switch dim
        case 2
            u0  = [u0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), t0, param )'];
            du0 = [du0; DATA.du0{k}( MESH.nodes(1,:), MESH.nodes(2,:), t0, param )'];
            
        case 3
            u0  = [u0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];
            du0 = [du0; DATA.du0{k}( MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];
    end
end
d2u0 = 0*du0;
u = u0;
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u0, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, 0);
end

TimeAdvance.Initialize( ROM.V'*u0(MESH.internal_dof), ROM.V'*du0(MESH.internal_dof), ROM.V'*d2u0(MESH.internal_dof) );
Coef_Mass = TimeAdvance.MassCoefficient( );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices         = %d \n', MESH.numVertices);
fprintf(' * Number of Nodes            = %d \n', MESH.numNodes);
fprintf(' * Number of Elements         = %d \n', MESH.numElem);
fprintf(' * Number of Reduced Elements = %d \n', ROM.Red_Mesh.numElem);
fprintf(' * %% Reduced Elements         = %2.2f \n', ROM.Red_Mesh.numElem / MESH.numElem * 100);
fprintf(' * Number of Dofs             = %d \n', size(ROM.V,2));
fprintf(' * Number of timesteps        = %d \n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

%% Preconditioner (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

SolidModel = CSM_Assembler( ROM.Red_Mesh, DATA, FE_SPACE );

%% Assemble mass matrix
fprintf('\n Assembling mass matrix... ');
t_assembly = tic;
M    =  DATA.Density * ROM.M;
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

LinSolver = LinearSolver( DATA.LinearSolver );

% Assemble matrix and right-hand side
fprintf('\n -- Assembling external Forces at t0... ');
t_assembly = tic;
F_ext_old_FE  = SolidModel.compute_external_forces(t0);
F_ext_old_FE = F_ext_old_FE( MESH.internal_dof );
F_ext_old = ROM.LeftProjection_ext * F_ext_old_FE( ROM.IDEIM_ext );
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

fprintf('\n -- Assembling internal Forces at t0... ');
t_assembly = tic;
F_in_old_FE  =  SolidModel.compute_internal_forces( u0 );
F_in_old_FE  = F_in_old_FE( MESH.internal_dof );
F_in_old     = ROM.LeftProjection_int * F_in_old_FE( ROM.IDEIM_in );
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly)

U_kN = ROM.V' * u(MESH.internal_dof,1);

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
    U_k            = u(:,end);
    U_k(MESH.Dirichlet_dof) = u_D;
    
    Csi = TimeAdvance.RhsContribute( );
    
    % Assemble matrix and right-hand side
    fprintf('\n -- Assembling external Forces... ');
    t_assembly = tic;
    F_ext_FE      = SolidModel.compute_external_forces( t );
    F_ext_FE      = F_ext_FE( MESH.internal_dof );
    F_ext         = ROM.LeftProjection_ext * F_ext_FE( ROM.IDEIM_ext );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    fprintf('\n -- Assembling internal Forces ... ');
    t_assembly = tic;
    F_in_FE      = SolidModel.compute_internal_forces( U_k );
    F_in_FE      = F_in_FE(MESH.internal_dof);
    F_in         = ROM.LeftProjection_int * F_in_FE( ROM.IDEIM_in );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    Residual  = Coef_Mass * M * U_kN +  ...
                  (1 - TimeAdvance.M_alpha_f) * F_in  + TimeAdvance.M_alpha_f * F_in_old ...
                - (1 - TimeAdvance.M_alpha_f) * F_ext + TimeAdvance.M_alpha_f * F_ext_old ...
                - M * Csi;
            
            
    fprintf('\n -- Assembling Jacobian matrix... ');
    t_assembly = tic;
    dF_in_FE   = SolidModel.compute_jacobian( U_k );
    dF_in_FE   =  CSM_ApplyBC(dF_in_FE, [], FE_SPACE, MESH, DATA, t, 1); % attenzione a carichi 
    dF_in      =  ROM.LeftProjection_int * ( dF_in_FE(ROM.IDEIM_in, : ) * ROM.V );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);        
            
    Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in;
    
%     % Apply boundary conditions
%     fprintf('\n -- Apply boundary conditions ... ');
%     t_assembly = tic;
%     [A_FE, b]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE, MESH, DATA, t, 1); % attenzione a carichi 
%     t_assembly = toc(t_assembly);
%     fprintf('done in %3.3f s\n', t_assembly);

    res0Norm = norm(Residual);
        
    fprintf('\n============ Start Newton Iterations ============\n\n');
    while (k <= maxIter && (incrNorm > tol || resRelNorm > tol) )
        
        % Solve
        fprintf('\n   -- Solve J x = -R ... ');
        Precon.Build( Jacobian );
        fprintf('\n        time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
        LinSolver.SetPreconditioner( Precon );
        dU_N                  = LinSolver.Solve( Jacobian, -Residual );
        dU(MESH.internal_dof) = ROM.V * dU_N;
        fprintf('\n        time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
        
        % update solution
        U_k        = U_k + dU;
        U_kN       = U_kN + dU_N;
        incrNorm   = norm(dU)/norm(U_k);
        
        % Assemble matrix and right-hand side
        fprintf('\n   -- Assembling internal forces... ');
        t_assembly = tic;
        F_in_FE      = SolidModel.compute_internal_forces( full ( U_k ) );
        F_in_FE      = F_in_FE(MESH.internal_dof);
        F_in         = ROM.LeftProjection_int * F_in_FE( ROM.IDEIM_in );
        dF_in_FE     = SolidModel.compute_jacobian( full ( U_k ) );
        dF_in_FE   =  CSM_ApplyBC(dF_in_FE, [], FE_SPACE, MESH, DATA, t, 1); % attenzione a carichi
        dF_in      =  ROM.LeftProjection_int * ( dF_in_FE(ROM.IDEIM_in, : ) * ROM.V );
        t_assembly = toc(t_assembly);
        fprintf('done in %3.3f s\n', t_assembly);
        
        Residual  = Coef_Mass * M * U_kN +  ...
                  (1 - TimeAdvance.M_alpha_f) * F_in  + TimeAdvance.M_alpha_f * F_in_old ...
                - (1 - TimeAdvance.M_alpha_f) * F_ext + TimeAdvance.M_alpha_f * F_ext_old ...
                - M * Csi;
            
        Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in;
        
        % Apply boundary conditions
%         fprintf('\n   -- Apply boundary conditions ... ');
%         t_assembly = tic;
%         [A, b]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE, MESH, DATA, t, 1);  % attenzione a carichi 
%         t_assembly = toc(t_assembly);
%         fprintf('done in %3.3f s\n', t_assembly);
        
        resRelNorm = norm(Residual) / res0Norm;
        
        fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k,full(incrNorm), full(norm(resRelNorm)));
        k = k + 1;
        
    end
    
    if sol_history
        u = [u U_k];
    else
        u = U_k;
    end
    
    %% Export to VTK
    if ~isempty(vtk_filename)
        CSM_export_solution(MESH.dim, U_k, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, k_t);
    end
    
    TimeAdvance.Update( U_kN );
    
    F_ext_old = F_ext;
    F_in_old  = F_in;
    
    iter_time = toc(iter_time);
    fprintf('\n-------------- Iteration time: %3.2f s -----------------',iter_time);
    
end

fprintf('\n************************************************************************* \n');

return
