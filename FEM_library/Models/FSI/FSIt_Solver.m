function [X, MESH, DATA] = FSIt_Solver(dim, meshFluid, meshSolid, fem_F, fem_S, data_file_F, data_file_S, param, vtk_filename)
%FSIT_SOLVER solves Fluid-Structure Interaction problems in 2D/3D

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 8
    param = [];
end

if nargin < 9
    vtk_filename = [];
end


%% Read problem parameters and BCs from data_file
DATA.Fluid   = CFD_read_DataFile(data_file_F, dim, param);
if nargin < 8
    DATA.Fluid.param = [];
else
    DATA.Fluid.param = param;
end

DATA.Solid   = CSM_read_DataFile(data_file_S, dim, param);
if nargin < 8
    DATA.Solid.param = [];
else
    DATA.Solid.param = param;
end

use_SUPG = false;
if isfield(DATA.Fluid, 'Stabilization')
    if strcmp( DATA.Fluid.Stabilization, 'SUPG' )
        use_SUPG = true;
    end
end

%% Set quad_order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

[ MESH.Fluid ] = buildMESH( dim, meshFluid.elements, meshFluid.vertices, ...
    meshFluid.boundaries, fem_F{1}, quad_order, DATA.Fluid, 'CFD', meshFluid.rings );
[ MESH.Solid ] = buildMESH( dim, meshSolid.elements, meshSolid.vertices, ...
    meshSolid.boundaries, fem_S, quad_order, DATA.Solid, 'CSM', meshSolid.rings );

MESH.dim  = dim;

%% Create and fill the FE_SPACE data structure
[ FE_SPACE_v ] = buildFESpace( MESH.Fluid, fem_F{1}, dim, quad_order );% fluid velocity
[ FE_SPACE_p ] = buildFESpace( MESH.Fluid, fem_F{2}, 1,   quad_order );% fluid pressure
[ FE_SPACE_s ] = buildFESpace( MESH.Solid, fem_S,    dim, quad_order );% solid displacement
[ FE_SPACE_g ] = buildFESpace( MESH.Fluid, fem_F{1}, dim, quad_order );% geometry displacement

MESH.Fluid.internal_dof_c{MESH.dim+1} = 1:FE_SPACE_p.numDof;

[MESH] = FSI_InterfaceMap(DATA, MESH);

for k = 1 : dim
    MESH.ndof_interface{k} = length(MESH.Interface_FSmap{k});
end

%% find interface indices in the internal numbering
tmp = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
for i = 1 : dim
    tmp([FE_SPACE_v.numDofScalar*(i-1)+MESH.Fluid.dof_interface{i}]) = 1;
end

MESH.Fluid.Gamma     = find(tmp(MESH.Fluid.internal_dof));
MESH.Fluid.II        = setdiff(1:length(MESH.Fluid.internal_dof),MESH.Fluid.Gamma);
MESH.Fluid.II_global = setdiff(MESH.Fluid.internal_dof,MESH.Fluid.internal_dof(MESH.Fluid.Gamma));

tmp = zeros(FE_SPACE_s.numDof,1);
for i = 1 : dim
    tmp([FE_SPACE_s.numDofScalar*(i-1)+MESH.Solid.dof_interface{i}]) = 1;
end
MESH.Solid.Gamma     = find(tmp(MESH.Solid.internal_dof));
MESH.Solid.II        = setdiff(1:length(MESH.Solid.internal_dof),MESH.Solid.Gamma);
MESH.Solid.II_global = setdiff(MESH.Solid.internal_dof,MESH.Solid.internal_dof(MESH.Solid.Gamma));

MESH.internal_dof    = [MESH.Fluid.II MESH.Fluid.Gamma' ...
            length(MESH.Fluid.internal_dof)+1:(length(MESH.Fluid.internal_dof)+length(MESH.Solid.II))];

Fluid_ReferenceNodes = MESH.Fluid.nodes(1:dim,:);

%% Gather Time Setting

BDF_orderF = DATA.Fluid.time.BDF_order;
t0        = DATA.Fluid.time.t0;
dt        = DATA.Fluid.time.dt;
tf        = DATA.Fluid.time.tf;
t         = DATA.Fluid.time.t0;

fprintf('\n **** FLUID PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n', MESH.Fluid.numVertices);
fprintf(' * Number of Elements  = %d \n', MESH.Fluid.numElem);
fprintf(' * Number of Nodes     = %d \n', MESH.Fluid.numNodes);
fprintf(' * Velocity DOFs       = %d \n', FE_SPACE_v.numDof);
fprintf(' * Pressure DOFs       = %d \n', FE_SPACE_p.numDof);
fprintf(' * BDF Order           = %d\n', BDF_orderF);


fprintf('\n **** SOLID PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n', MESH.Solid.numVertices);
fprintf(' * Number of Elements  = %d \n', MESH.Solid.numElem);
fprintf(' * Number of Nodes     = %d \n', MESH.Solid.numNodes);
fprintf(' * Displacement DOFs   = %d \n', FE_SPACE_s.numDof);

fprintf('\n * FSI interface DOFs  =  %d\n', sum(cell2mat(MESH.ndof_interface)));
fprintf(' * Linear System dim   =  %d\n', length(MESH.internal_dof));
fprintf(' * Number of timesteps =  %d\n', (tf-t0)/dt);
fprintf('-------------------------------------------\n');

X_n  = zeros(length(MESH.Fluid.II) + length(MESH.Fluid.Gamma) + length(MESH.Solid.II),1);
X_nk = X_n;

k_t       = 0;

%% Initalize Fluid Time Advance
TimeAdvanceF = BDF_TimeAdvance( BDF_orderF );

v0  = [];
for k = 1 : FE_SPACE_v.numComponents
    switch dim
        case 2
            v0  = [v0; DATA.Fluid.u0{k}(  MESH.Fluid.nodes(1,:), MESH.Fluid.nodes(2,:), t0, param )'];
            
        case 3
            v0  = [v0; DATA.Fluid.u0{k}(  MESH.Fluid.nodes(1,:), MESH.Fluid.nodes(2,:), MESH.Fluid.nodes(3,:), t0, param )'];
    end
end
u = [v0; zeros(FE_SPACE_p.numDof,1)];
if ~isempty(vtk_filename)
    CFD_export_solution(MESH.dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), ...
        MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], 0);
end

TimeAdvanceF.Initialize( v0 );

%% Initalize Solid Time Advance
TimeAdvanceS = Newmark_TimeAdvance( DATA.Solid.time.beta, DATA.Solid.time.gamma, dt );

u0  = [];
du0 = [];
for k = 1 : FE_SPACE_s.numComponents
    switch dim
        case 2
            u0  = [u0; DATA.Solid.u0{k}(  MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), t0, param )'];
            du0 = [du0; DATA.Solid.du0{k}( MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), t0, param )'];
            
        case 3
            u0  = [u0; DATA.Solid.u0{k}(  MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), MESH.Solid.nodes(3,:), t0, param )'];
            du0 = [du0; DATA.Solid.du0{k}( MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), MESH.Solid.nodes(3,:), t0, param )'];
    end
end
d2u0 = 0*du0;
u = u0;
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u0, MESH.Solid.vertices, MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename,'Solid'], 0);
end

TimeAdvanceS.Initialize( u0, du0, d2u0 );
Coef_MassS = TimeAdvanceS.MassCoefficient( );


%% Initalize Geometry Time Advance
ALE_velocity = zeros(MESH.Fluid.numNodes*dim, 1);
%TimeAdvanceG = BDF_TimeAdvance( BDF_orderF );
%TimeAdvanceG.Initialize( ALE_velocity );

d_Fn      = zeros(MESH.Fluid.numNodes*dim, 1);

%% Generate Domain Decomposition (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Fluid.Preconditioner.type, DATA.Fluid);

% if isfield(DATA.Preconditioner, 'type') && strcmp( DATA.Preconditioner.type, 'AdditiveSchwarz')
%     R      = FSI_overlapping_DD(MESH, FE_SPACE_v, FE_SPACE_p, DATA.Preconditioner.num_subdomains,  DATA.Preconditioner.overlap_level);
%     Precon.SetRestrictions( R );
% end

%% assemble Harmonic Extension matrix
fprintf('\n Assemble Harmonic Extension matrix and store LU factors ... '); 
time2          = tic; 

DATA.Geometry                  = DATA.Solid;
DATA.Geometry.Material_Model   = 'SEMMT';
DATA.Geometry.Stiffening_power = 0.8;

MeshMotionAssembler = CSM_Assembler( MESH.Fluid, DATA.Geometry, FE_SPACE_g );

Harmonic_Ext_Matrix.HE = MeshMotionAssembler.compute_jacobian( zeros(FE_SPACE_g.numDof, 1) );

internal_dofs_HE = [];
for k = 1 : dim
    internal_dofs_HE       = [ internal_dofs_HE (k-1)*FE_SPACE_v.numDofScalar+setdiff(1:FE_SPACE_v.numDofScalar, ...
        [MESH.Fluid.dof_interface{k}; MESH.ALE_dirichlet{k}])];
end
Harmonic_Ext_Matrix.internal_dofs = internal_dofs_HE;
[Harmonic_Ext_Matrix.L , Harmonic_Ext_Matrix.U,...
    Harmonic_Ext_Matrix.perm , q ]   = lu(Harmonic_Ext_Matrix.HE(internal_dofs_HE,internal_dofs_HE), 'vector');
Harmonic_Ext_Matrix.invp             = 0*q ;
Harmonic_Ext_Matrix.invp(q)          = 1:length(q);

time2          = toc(time2);
fprintf('%f s \n',time2);
fprintf('-------------------------------------------\n');

%% Coupling matrices
Z_FS  = sparse(length(MESH.Fluid.II), length(MESH.Solid.II));
Z_SF  = sparse(length(MESH.Solid.II), length(MESH.Fluid.II));

IdGamma_SF = [];

for k = 1 : dim
    if ~isempty(MESH.ndof_interface{k})
        IdGamma_SF_tmp = sparse(MESH.ndof_interface{k}, MESH.ndof_interface{k});
        IdGamma_SF_tmp(MESH.Interface_FSmap{k}, :) = speye(MESH.ndof_interface{k},MESH.ndof_interface{k});
    else
        IdGamma_SF_tmp = [];
    end
    IdGamma_SF = blkdiag(IdGamma_SF, IdGamma_SF_tmp);
end

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

%% Create Solid Assembler Object
SolidModel = CSM_Assembler( MESH.Solid, DATA.Solid, FE_SPACE_s );

% Assemble mass matrix
fprintf('\n Assembling S-mass matrix... ');
t_assembly = tic;
M_s    =  SolidModel.compute_mass();
M_s    =  M_s * DATA.Solid.Density;
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

if strcmp( DATA.Solid.Material_Model, 'Linear' )
    A_s = SolidModel.compute_jacobian( zeros(FE_SPACE_s.numDof, 1) );
end

%% Initialize Linear Solver
LinSolver = LinearSolver( DATA.Fluid.LinearSolver );

tol        = DATA.Solid.NonLinearSolver.tol;
maxIter    = DATA.Solid.NonLinearSolver.maxit;


%% PreProcessing for Drag and Lift Computation
compute_AerodynamicForces = 0 ;
if isfield(DATA.Fluid, 'Output') && isfield(DATA.Fluid.Output, 'DragLift')
    if DATA.Fluid.Output.DragLift.computeDragLift == 1
        compute_AerodynamicForces = true;
    end
end

if compute_AerodynamicForces
    AeroF_x(k_t+1)  = 0;
    AeroF_y(k_t+1)  = 0;
    AeroF_z(k_t+1)  = 0;
    dofs_drag    = [];
    
    for j = 1 : length(DATA.Fluid.Output.DragLift.flag)
        Dirichlet_side         = find(MESH.Fluid.boundaries(MESH.Fluid.bc_flag_row,:) == DATA.Fluid.Output.DragLift.flag(j));
        Dirichlet_side         = unique(Dirichlet_side);
        Dirichlet_dof          = MESH.Fluid.boundaries(1:MESH.Fluid.numBoundaryDof,Dirichlet_side);
        dofs_drag              = [dofs_drag; Dirichlet_dof(:)];
    end
    dofs_drag = unique(dofs_drag);
    
    fileDragLift = fopen(DATA.Fluid.Output.DragLift.filename, 'w+');
    fprintf(fileDragLift, 'Time          F_x          F_y          F_z');
    fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
end


%% Time Loop
fprintf('\n **** Starting temporal loop ****\n');
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n',t0,t,tf);
    
    v_BDF = TimeAdvanceF.RhsContribute( );
    u_BDF = [v_BDF; zeros(FE_SPACE_p.numDof,1)];
    alpha = TimeAdvanceF.GetCoefficientDerivative();
    
    switch DATA.Fluid.time.nonlinearity
        
        case 'semi-implicit'
            
            %% Update Fluid Linear Matrices
            FluidModel = CFD_Assembler( MESH.Fluid, DATA.Fluid, FE_SPACE_v, FE_SPACE_p );

            fprintf('\n   -- Fluid_Assembling Stokes terms... ');
            t_assembly = tic;
            [A_Stokes] = FluidModel.compute_Stokes_matrix();
            t_assembly = toc(t_assembly);
            fprintf('done in %3.3f s\n', t_assembly);
            
            fprintf('\n   -- Fluid_Assembling mass matrix... ');
            t_assembly = tic;
            Mv = FluidModel.compute_mass_velocity();
            Mp = FluidModel.compute_mass_pressure();
            M  = blkdiag(DATA.Fluid.density * Mv, 0*Mp);
            t_assembly = toc(t_assembly);
            fprintf('done in %3.3f s\n', t_assembly);

            v_extrapolated = TimeAdvanceF.Extrapolate();
            
            fprintf('\n   -- Fluid_Assembling Convective Term... ');
            t_assembly = tic;
            [C1] = FluidModel.compute_convective_Oseen_matrix( v_extrapolated - ALE_velocity);
            t_assembly = toc(t_assembly);
            fprintf('done in %3.3f s\n', t_assembly);
            
            F_NS = 1/dt * M * u_BDF;
            C_NS = alpha/dt * M + A_Stokes + C1;
            
            if use_SUPG
                fprintf('\n   -- Fluid_Assembling SUPG Terms ... ');
                t_assembly = tic;
                [A_SUPG, F_SUPG] = FluidModel.compute_SUPG_semiimplicit( v_extrapolated - ALE_velocity, v_BDF, dt, alpha);
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                
                C_NS             = C_NS + A_SUPG;
                F_NS             = F_NS - F_SUPG;
            end
            
            % Apply Fluid boundary conditions
            fprintf('\n   -- Fluid_Apply boundary conditions ... ');
            t_assembly = tic;
            [C_NS_in, F_NS_in, v_D]   =  CFD_ApplyBC(C_NS, F_NS, FE_SPACE_v, FE_SPACE_p, MESH.Fluid, DATA.Fluid, t);
            t_assembly = toc(t_assembly);
            fprintf('done in %3.3f s\n', t_assembly);
            
            switch DATA.Solid.Material_Model
                
                case 'Linear'                    
                    
                    Csi           = TimeAdvanceS.RhsContribute( );
                    F_ext         = SolidModel.compute_volumetric_forces( t );
                    
                    % Solid matrix and rhs Update
                    C_STR  = Coef_MassS * M_s +  A_s;
                    F_S    = F_ext + M_s * Csi;
                    [C_STR, F_S, DisplacementDir_np1] = CSM_ApplyBC(C_STR, F_S, FE_SPACE_s, MESH.Solid, DATA.Solid, t);
                    
                    % Get displacment d^alpha on the interface
                    %d_alpha = dDisplacement_n - dt * DATA.Solid.time.gamma * Csi + dt* (1-DATA.Solid.time.gamma)*d2Displacement_n;
                    d_alpha = TimeAdvanceS.M_dU - dt * DATA.Solid.time.gamma * Csi + dt* (1-DATA.Solid.time.gamma)*TimeAdvanceS.M_d2U;

                    F_L     = d_alpha(MESH.Solid.internal_dof(MESH.Solid.Gamma));
                    alpha   = DATA.Solid.time.gamma / (dt * DATA.Solid.time.beta);
                    
                    % Monolothic System
                    S_GG     = IdGamma_SF * (C_STR(MESH.Solid.Gamma,MESH.Solid.Gamma) * IdGamma_FS);
                    S_GI     = IdGamma_SF *  C_STR(MESH.Solid.Gamma,MESH.Solid.II);
                    
                    FSI_M    = [C_NS_in(MESH.Fluid.II,MESH.Fluid.II)                      C_NS_in(MESH.Fluid.II,MESH.Fluid.Gamma)            Z_FS ;...
                        C_NS_in(MESH.Fluid.Gamma,MESH.Fluid.II)      C_NS_in(MESH.Fluid.Gamma,MESH.Fluid.Gamma)+1/alpha*S_GG                 S_GI  ;...
                        Z_SF                                                    1/alpha*C_STR(MESH.Solid.II,MESH.Solid.Gamma)*IdGamma_FS     C_STR(MESH.Solid.II,MESH.Solid.II)];
                    
                    FSI_F    = [F_NS_in(MESH.Fluid.II); ...
                        F_NS_in(MESH.Fluid.Gamma)+IdGamma_SF*F_S(MESH.Solid.Gamma)+1/alpha*(IdGamma_SF*(C_STR(MESH.Solid.Gamma,MESH.Solid.Gamma)*F_L)); ...
                        F_S(MESH.Solid.II)+1/alpha*C_STR(MESH.Solid.II,MESH.Solid.Gamma)*F_L];
                                        
                    % Solve
                    fprintf('\n -- Solve A x = b ... ');
                    Precon.Build( FSI_M );
                    fprintf('\n      time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
                    LinSolver.SetPreconditioner( Precon );
                    X_nk(MESH.internal_dof) = LinSolver.Solve( FSI_M, FSI_F, X_n(MESH.internal_dof) );
                    fprintf('\n      time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
                    
                otherwise
                    error('to be coded')
                    
            end
            
        
    end
    
    norm_n   = norm(X_nk - X_n) / norm(X_n);
    X_n      = X_nk;

    %% Export solid displacement on reference mesh
    Displacement_np1                                   = zeros(FE_SPACE_s.numDof,1);
    Displacement_np1(MESH.Solid.II_global)             = X_n(1+length(MESH.Fluid.internal_dof):end);
    Displacement_np1(MESH.Solid.Dirichlet_dof)         = DisplacementDir_np1;
    Displacement_np1(MESH.Solid.internal_dof(MESH.Solid.Gamma))   = 1/alpha * ( IdGamma_FS*X_n(MESH.Fluid.Gamma) -  F_L);
    
    TimeAdvanceS.Update( Displacement_np1 );

    % Export to VTK
    if ~isempty(vtk_filename)
        CSM_export_solution(MESH.dim, Displacement_np1, MESH.Solid.vertices, ...
            MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename, 'Solid'], k_t);
    end
    %d2Displacement_np1      = 1 / ( TIMEs.beta * dt^2) * Displacement_np1 - Csi;
    %dDisplacement_np1       = dDisplacement_n + dt * (TIMEs.gamma * d2Displacement_np1 + (1 - TIMEs.gamma) * d2Displacement_n) ;
    %Displacement_n          = Displacement_np1;
    %dDisplacement_n         = dDisplacement_np1;
    %d2Displacement_n        = d2Displacement_np1;
    
    %% Deform Fluid mesh by Solid-Extension
    
    d_F = FSI_harmonicExtension(MESH, Displacement_np1, Harmonic_Ext_Matrix);
    Fluid_def_nodes    = Fluid_ReferenceNodes + d_F;
    Fluid_def_vertices = Fluid_def_nodes(1:dim, 1:MESH.Fluid.numVertices);
        
    d_F      = reshape(d_F',dim*MESH.Fluid.numNodes,1);
    
    %% Compute Fluid mesh velocity: w = 1/dt * ( d_f^(n+1) - d_f^n )
    ALE_velocity  =  1/dt * ( d_F - d_Fn );
    d_Fn          =  d_F;
    
    %% Update Fluid MESH
    MESH.Fluid.vertices = Fluid_def_vertices;
    MESH.Fluid.nodes    = Fluid_def_nodes;
    [MESH.Fluid.jac, MESH.Fluid.invjac, MESH.Fluid.h] = geotrasf(dim, MESH.Fluid.vertices, MESH.Fluid.elements);
    
    %% Export Fluid velocity and pressure on deformed mesh
    u(MESH.Fluid.internal_dof)  = X_n(1:length(MESH.Fluid.internal_dof));
    u(MESH.Fluid.Dirichlet_dof) = v_D;
        
    % Export to VTK
    if ~isempty(vtk_filename)
        CFD_export_solution(dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), ...
            MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], k_t);
    end
    
    %% Update BDF
    TimeAdvanceF.Append( u(1:FE_SPACE_v.numDof) );
    
    %% Compute_DragLift
    if compute_AerodynamicForces
        
        if strcmp(DATA.Fluid.time.nonlinearity,'implicit')
            %C_NS = 0*Jacobian;
            %F_NS = -Residual;
            error('Areo Forces computation to be coded')
        end
        Z              = zeros(FE_SPACE_v.numDofScalar,1);
        Z(dofs_drag)   = 1;
        
        W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
        W(1:FE_SPACE_v.numDofScalar)        = Z;
        AeroF_x(k_t+1) = DATA.Fluid.Output.DragLift.factor*(W'*(-C_NS*u + F_NS));
        
        W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
        W(FE_SPACE_v.numDofScalar+[1:FE_SPACE_v.numDofScalar])  = Z;
        AeroF_y(k_t+1)  = DATA.Fluid.Output.DragLift.factor*(W'*(-C_NS*u  + F_NS));
        
        if MESH.dim == 3
            W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
            W(2*FE_SPACE_v.numDofScalar+[1:FE_SPACE_v.numDofScalar])  = Z;
            AeroF_z(k_t+1)  = DATA.Fluid.Output.DragLift.factor*(W'*(-C_NS*u  + F_NS));
        else
            AeroF_z(k_t+1) = 0.0;
        end
        
        fprintf('\n *** F_x = %e, F_y = %e, F_z = %e *** \n',  AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
        fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
    end

    iter_time = toc(iter_time);
    fprintf('\nnorm(U^(n+1) - U^(n))/norm(U^(n)) = %2.3e -- Iteration time: %3.2f s \n',full(norm_n),iter_time);
    %fprintf('\n-------------- Iteration time: %3.2f s -----------------',iter_time);
    
    %X(:,k_t+1) = [u; Displacement_np1];
    X = [u; Displacement_np1];
    
end

fprintf('\n************************************************************************* \n');

return
