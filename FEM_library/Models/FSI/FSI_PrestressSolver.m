function [R_P, J_P] = FSI_PrestressSolver(dim, meshFluid, meshSolid, fem_F, fem, data_file_F, data_file_S, param, vtk_filename)
%FSI_PRESTRESSSOLVER given a fluid load on the FS interface, solves modified 
%nonlinear elastodynamics equations to find the corresponding prestress
%tensor, following the procedure detailed in:
%
%   Hsu, M. C., & Bazilevs, Y. (2011). Blood vessel tissue prestress modeling
%   for vascular fluid?structure interaction simulation. Finite Elements in 
%   Analysis and Design, 47(6), 593-599.

%
%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

%% Read problem parameters and BCs from data_file
DATA.Solid   = CSM_read_DataFile(data_file_S, dim, param);
DATA.Fluid   = CFD_read_DataFile(data_file_F, dim, param);

DATA.Solid.param = param;
DATA.Fluid.param = param;

%% Set quad_order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

%% Create and fill the MESH data structure
[ MESH.Solid ] = buildMESH( dim, meshSolid.elements, meshSolid.vertices, ...
    meshSolid.boundaries, fem, quad_order, DATA.Solid, 'CSM', meshSolid.rings );

[ MESH.Fluid ] = buildMESH( dim, meshFluid.elements, meshFluid.vertices, ...
    meshFluid.boundaries, fem_F{1}, quad_order, DATA.Fluid, 'CFD', meshFluid.rings );

MESH.dim  = dim;

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH.Solid, fem, dim, quad_order );

[ FE_SPACE_v ] = buildFESpace( MESH.Fluid, fem_F{1}, dim, quad_order );% fluid velocity
[ FE_SPACE_p ] = buildFESpace( MESH.Fluid, fem_F{2}, 1,   quad_order );% fluid pressure
MESH.Fluid.internal_dof_c{MESH.dim+1} = 1:FE_SPACE_p.numDof;

%% Generates mappings from solid to fluid interface dofs and viceversa
fprintf('\n Build FS interface Maps... ');
t_assembly = tic;
[MESH] = FSI_InterfaceMap(DATA, MESH);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

for k = 1 : dim
    MESH.ndof_interface{k} = length(MESH.Interface_FSmap{k});
end

load(DATA.Fluid.Output.ExportFinalLoad.filename); % load FluidLoad

% restrict FluidLoad to FS interface
Fluid_DofsLoad = [];
for i = 1 : dim
    Fluid_DofsLoad = [Fluid_DofsLoad; FE_SPACE_v.numDofScalar*(i-1)+MESH.Fluid.dof_interface{i}];
end

Solid_DofsLoad = [];
for i = 1 : dim
    Solid_DofsLoad = [Solid_DofsLoad; FE_SPACE.numDofScalar*(i-1)+MESH.Solid.dof_interface{i}];
end

% transfer matrix from fluid to solid
IdGamma_FS = [];
for k = 1 : dim
    if ~isempty(MESH.ndof_interface{k})
        IdGamma_FS_tmp = sparse(MESH.ndof_interface{k}, MESH.ndof_interface{k});
        IdGamma_FS_tmp(MESH.Interface_SFmap{k}, :) = speye(MESH.ndof_interface{k},MESH.ndof_interface{k});
    else
        IdGamma_FS_tmp = [];
    end
    IdGamma_FS = blkdiag(IdGamma_FS, IdGamma_FS_tmp);
end

StructureLoad                 = zeros(FE_SPACE.numDof, 1);
StructureLoad(Solid_DofsLoad) = IdGamma_FS * FluidLoad(Fluid_DofsLoad);
F_extFluid                    = @(t) DATA.Solid.LoadTimeProfile(t) * StructureLoad;

% Test Exporting
StructureLoadF                 = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof, 1);
StructureLoadF(Fluid_DofsLoad) = FluidLoad(Fluid_DofsLoad);

CSM_export_solution(MESH.Solid.dim, StructureLoad, MESH.Solid.vertices, MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename, 'SLoad'], 0);
CFD_export_solution(MESH.dim, StructureLoadF(1:FE_SPACE_v.numDof), StructureLoadF(1+FE_SPACE_v.numDof:end), MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename, 'FLoad'], 0);

 
%% Gather Time Setting
t0        = DATA.Solid.time.t0;
dt        = DATA.Solid.time.dt;
tf        = DATA.Solid.time.tf;
t         = DATA.Solid.time.t0;
k_t       = 0;

TimeAdvance = GeneralizedAlpha_TimeAdvance( DATA.Solid.time.beta, DATA.Solid.time.gamma, DATA.Solid.time.alpha_m, DATA.Solid.time.alpha_f, dt );

u0  = [];
du0 = [];
for k = 1 : FE_SPACE.numComponents
    switch dim
        case 2
            u0  = [u0; DATA.Solid.u0{k}(  MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), t0, param )'];
            du0 = [du0; DATA.Solid.du0{k}( MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), t0, param )'];
            
        case 3
            u0  = [u0; DATA.Solid.u0{k}(  MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), MESH.Solid.nodes(3,:), t0, param )'];
            du0 = [du0; DATA.Solid.du0{k}( MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), MESH.Solid.nodes(3,:), t0, param )'];
    end
end

u = u0;
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.Solid.dim, u0, MESH.Solid.vertices, MESH.Solid.elements, MESH.Solid.numNodes, vtk_filename, 0);
end

Coef_Mass = TimeAdvance.MassCoefficient( );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.Solid.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.Solid.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.Solid.numNodes);
fprintf(' * Number of Dofs      = %d \n',length(MESH.Solid.internal_dof));
fprintf(' * Number of timesteps =  %d\n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

%% Generate Domain Decomposition (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Solid.Preconditioner.type, DATA.Solid);

if isfield(DATA.Solid.Preconditioner, 'type') && strcmp( DATA.Solid.Preconditioner.type, 'AdditiveSchwarz')
    R      = CSM_overlapping_DD(MESH.Solid, DATA.Solid.Preconditioner.num_subdomains,  DATA.Solid.Preconditioner.overlap_level);
    Precon.SetRestrictions( R );
end

SolidModel = CSM_Assembler( MESH.Solid, DATA.Solid, FE_SPACE );

%% Assemble mass matrix
fprintf('\n Assembling mass matrix... ');
t_assembly = tic;
M    =  SolidModel.compute_mass();
M    =  M * DATA.Solid.Density;
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

% Assemble Robin BC (if it's the case)
A_robin = SolidModel.assemble_ElasticRobinBC();

%% Zero Initial Acceleration
d2u0 = 0 * u0;

TimeAdvance.Initialize( u0, du0, d2u0 );

LinSolver = LinearSolver( DATA.Solid.LinearSolver );

U_n = u0;

%% Initialize Prestress
[R_P, J_P, S_np1] = SolidModel.compute_prestress(U_n, []);

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );
    fprintf('\n==========  Load percentage = %2.3f\n', DATA.Solid.LoadTimeProfile(t) * 100 );
    
    % Newton Method
    tol        = DATA.Solid.NonLinearSolver.tol;
    resRelNorm = tol + 1;
    incrNorm   = tol + 1;
    maxIter    = DATA.Solid.NonLinearSolver.maxit;
    k          = 1;
    
    [~, ~, u_D]   =  CSM_ApplyBC([], [], FE_SPACE, MESH.Solid, DATA.Solid, t);
    dU             = zeros(MESH.Solid.numNodes*MESH.Solid.dim,1);
    U_k            = u(:,end);
    U_k(MESH.Solid.Dirichlet_dof) = u_D;
    
    Csi = TimeAdvance.RhsContribute( );
    
    % Assemble matrix and right-hand side
    fprintf('\n -- Assembling external Forces... ');
    t_assembly = tic;
    F_ext      = SolidModel.compute_volumetric_forces( (1 - TimeAdvance.M_alpha_f) * t + TimeAdvance.M_alpha_f * (t-dt) );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    fprintf('\n -- Assembling internal Forces ... ');
    t_assembly = tic;
    F_in      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    Residual  = Coef_Mass * M * U_k + F_in - F_ext - M * Csi ...
                + A_robin * ((1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n) ...
                + F_extFluid( (1 - TimeAdvance.M_alpha_f) * t + TimeAdvance.M_alpha_f * (t-dt) ) ...
                + R_P + J_P * ((1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n);
            
    fprintf('\n -- Assembling Jacobian matrix... ');
    t_assembly = tic;
    dF_in     = SolidModel.compute_jacobian( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);        
            
    Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in ...
                + A_robin * (1 - TimeAdvance.M_alpha_f) ...
                + (1 - TimeAdvance.M_alpha_f) * J_P;
    
    % Apply boundary conditions
    fprintf('\n -- Apply boundary conditions ... ');
    t_assembly = tic;
    [A, b]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE, MESH.Solid, DATA.Solid, t, 1);
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
        dU(MESH.Solid.internal_dof) = LinSolver.Solve( A, b );
        fprintf('\n        time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
        
        % update solution
        U_k        = U_k + dU;
        incrNorm   = norm(dU)/norm(U_k);
        
        % Assemble matrix and right-hand side
        fprintf('\n   -- Assembling internal forces... ');
        t_assembly = tic;
        F_in       = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
        dF_in      = SolidModel.compute_jacobian( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
        t_assembly = toc(t_assembly);
        fprintf('done in %3.3f s\n', t_assembly);
        
        Residual  = Coef_Mass * M * U_k + F_in - F_ext - M * Csi ...
                    + A_robin * ((1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n) ...
                    + F_extFluid( (1 - TimeAdvance.M_alpha_f) * t + TimeAdvance.M_alpha_f * (t-dt) ) ...
                    + R_P + J_P * ((1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n);
            
        Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in ...
                    + A_robin * (1 - TimeAdvance.M_alpha_f) ...
                    + (1 - TimeAdvance.M_alpha_f) * J_P;
        
        % Apply boundary conditions
        fprintf('\n   -- Apply boundary conditions ... ');
        t_assembly = tic;
        [A, b]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE, MESH.Solid, DATA.Solid, t, 1);
        t_assembly = toc(t_assembly);
        fprintf('done in %3.3f s\n', t_assembly);
        
        resRelNorm = norm(b) / res0Norm;
        
        fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k,full(incrNorm), full(norm(resRelNorm)));
        k = k + 1;
        
    end
    fprintf('\n -- Norm(U_np1 - U_n) / Norm( U_n ) = %1.2e, Norm( U_np1 ) = %1.2e \n', norm(U_k - u) / norm(u), norm(U_k));
        
    %% Export to VTK
    if ~isempty(vtk_filename)
        CSM_export_solution(MESH.Solid.dim, U_k, MESH.Solid.vertices, MESH.Solid.elements, MESH.Solid.numNodes, vtk_filename, k_t);
    end
    
    %% Compute Von Mises Stress
    if DATA.Solid.Output.ComputeVonMisesStress
        fprintf('\n   -- Compute Element Stresses... ');
        t_assembly = tic;
        [~, Sigma]  =  SolidModel.compute_stress(U_k);
        t_assembly = toc(t_assembly);
        fprintf('done in %3.3f s\n', t_assembly);
        
        if MESH.Solid.dim == 2
            Sigma_VM = sqrt(  Sigma(:,1).^2 + Sigma(:,4).^2 - Sigma(:,1) .* Sigma(:,4) + 3 * Sigma(:,2).^2 );
        elseif MESH.Solid.dim == 3
            Sigma_VM = sqrt( 0.5 * ( (Sigma(:,1) - Sigma(:,5)).^2 + (Sigma(:,5) - Sigma(:,9)).^2 + (Sigma(:,9) - Sigma(:,1)).^2 + 6 * ( Sigma(:,2).^2 + Sigma(:,6).^2 + Sigma(:,7).^2 ) ) );
        end
        CSM_export_VonMisesStress(MESH.Solid.dim, Sigma_VM, MESH.Solid.vertices, MESH.Solid.elements, [vtk_filename, 'VMstress_'], k_t);
    end
    
    u   = U_k;
    TimeAdvance.Update( U_k );
    
    U_n = U_k;
    
    % Update Prestress
    fprintf('\n   -- Update Prestress ... ');
    t_assembly = tic;
    [R_P, J_P, S_np1] = SolidModel.compute_prestress(U_n, S_np1);
    t_assembly = toc(t_assembly);
    fprintf('done in %3.3f s\n', t_assembly);
    
    iter_time = toc(iter_time);
    fprintf('\n-------------- Iteration time: %3.2f s -----------------',iter_time);
    
end

fprintf('\n************************************************************************* \n');

return
