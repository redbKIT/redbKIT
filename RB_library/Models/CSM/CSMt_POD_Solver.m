function [u, MassMatrix, FE_SPACE, MESH, DATA] = CSMt_POD_Solver(dim, elements, vertices, boundaries, fem, data_file, ...
    param, vtk_filename, sol_history, Training_Options, V_POD)
%CSMT_POD_SOLVER Dynamic Structural POD-Galerkin Solver
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

u = u0;
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u0, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, 0);
end

Coef_Mass = TimeAdvance.MassCoefficient( );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',size(V_POD, 2));
fprintf(' * Number of timesteps =  %d\n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

if export_h5
    Fint_h5  = HDF5_DenseMultiCVector(Training_Options.h5_filename, ...
        Training_Options.InternalForces.h5_section, length(MESH.internal_dof));
    Fext_h5  = HDF5_DenseMultiCVector(Training_Options.h5_filename, ...
        Training_Options.ExternalForces.h5_section, length(MESH.internal_dof));
end

%% Preconditioner (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Preconditioner.type, DATA);

SolidModel = CSM_Assembler( MESH, DATA, FE_SPACE );

%% Assemble mass matrix
fprintf('\n Assembling mass matrix... ');
t_assembly = tic;
MassMatrix = SolidModel.compute_mass();
M_h        = MassMatrix * DATA.Density;
M_h_in     = M_h(MESH.internal_dof, MESH.internal_dof);
M          = V_POD' * (M_h_in * V_POD);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

LinSolver = LinearSolver( DATA.LinearSolver );

%% Initial Acceleration
fprintf('\n -- Assembling external Forces at t0... ');
t_assembly = tic;
F_ext_0_FE = SolidModel.compute_volumetric_forces(t0);
F_ext_0    = V_POD' * F_ext_0_FE( MESH.internal_dof );
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly);

if export_h5
    Fext_h5.append( F_ext_0_FE(MESH.internal_dof) );
end

fprintf('\n -- Assembling internal Forces at t0... ');
t_assembly = tic;
F_in_0_FE  = SolidModel.compute_internal_forces( u0 );
F_in_0     = V_POD' * F_in_0_FE( MESH.internal_dof );
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s\n', t_assembly)

if export_h5
    Fint_h5.append( F_in_0_FE(MESH.internal_dof) );
end

d2u0_N = M \ (F_ext_0 - F_in_0);

TimeAdvance.Initialize( V_POD'*u0(MESH.internal_dof), V_POD'*du0(MESH.internal_dof), d2u0_N);%ROM.V'*d2u0(MESH.internal_dof) );

U_kN = V_POD'*u0(MESH.internal_dof);
U_n  = u0;

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
    
    [~, ~, u_D]   =  CSM_ApplyEssentialBC([], [], MESH, DATA, t);
    dU             = zeros(MESH.numNodes*MESH.dim,1);
    U_k            = u(:,end);
    U_k(MESH.Dirichlet_dof) = u_D;
    
    Csi = TimeAdvance.RhsContribute( );
    
    % Assemble matrix and right-hand side
    fprintf('\n -- Assembling external Forces... ');
    t_assembly = tic;
    F_ext_FE      = SolidModel.compute_external_forces( (1 - TimeAdvance.M_alpha_f) * t + TimeAdvance.M_alpha_f * (t-dt) );
    F_ext         = V_POD' * F_ext_FE( MESH.internal_dof );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    if export_h5
        Fext_h5.append( F_ext_FE(MESH.internal_dof) );
    end
    
    fprintf('\n -- Assembling internal Forces ... ');
    t_assembly = tic;
    F_in_FE      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
    F_in         = V_POD' * F_in_FE( MESH.internal_dof );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    if export_h5
        Fint_h5.append( F_in_FE(MESH.internal_dof) );
    end
    
    Residual  = Coef_Mass * M * U_kN + F_in - F_ext - M * Csi;
            
    fprintf('\n -- Assembling Jacobian matrix... ');
    t_assembly = tic;
    dF_in_FE   = SolidModel.compute_jacobian( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
    dF_in_FE   = CSM_ApplyEssentialBC(dF_in_FE, [], MESH, DATA, t, 1);
    dF_in      = V_POD' * ( dF_in_FE * V_POD );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);        
            
    Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in;
    
    res0Norm  = norm(Residual);
        
    fprintf('\n============ Start Newton Iterations ============\n\n');
    while (k <= maxIter && (incrNorm > tol || resRelNorm > tol))
        
        % Solve
        fprintf('\n   -- Solve J x = -R ... ');
        Precon.Build( Jacobian );
        fprintf('\n        time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
        LinSolver.SetPreconditioner( Precon );
        dU_N                  = LinSolver.Solve( Jacobian, -Residual );
        dU(MESH.internal_dof) = V_POD * dU_N;
        fprintf('\n        time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
      
        
        % update solution
        U_k        = U_k + dU;
        U_kN       = U_kN + dU_N;
        incrNorm   = norm(dU)/norm(U_k);
        
        % Assemble matrix and right-hand side
        fprintf('\n   -- Assembling internal forces... ');
        t_assembly = tic;
        F_in_FE      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
        F_in         = V_POD' * F_in_FE(MESH.internal_dof);
        dF_in_FE     = SolidModel.compute_jacobian( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
        dF_in_FE     = CSM_ApplyEssentialBC(dF_in_FE, [], MESH, DATA, t, 1);
        dF_in        = V_POD' * ( dF_in_FE * V_POD );
        t_assembly = toc(t_assembly);
        fprintf('done in %3.3f s\n', t_assembly);
        
        if export_h5
            Fint_h5.append( F_in_FE(MESH.internal_dof) );
        end
        
        Residual  = Coef_Mass * M * U_kN + F_in - F_ext - M * Csi;
    
        Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in;
        
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
        
    U_n = U_k;
    
    iter_time = toc(iter_time);
    fprintf('\n-------------- Iteration time: %3.2f s -----------------',iter_time);
    
end

fprintf('\n************************************************************************* \n');

return
