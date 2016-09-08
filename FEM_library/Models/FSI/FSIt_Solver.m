function [X, MESH, DATA] = FSIt_Solver(dim, meshFluid, meshSolid, fem_F, fem_S, data_file_F, data_file_S, param, vtk_filename)
%FSIT_SOLVER solves Fluid-Structure Interaction problems in 2D/3D
%
%   The solver implements the following numerical approximations of the FSI
%   problem:
%
%   ======================================================================
%   ** OPT1: monolithic GCE with semi-implicit fluid and 
%            linear/nonlinear structure **
%  
%   This options is enabled by setting 
%               data.time.nonlinearity  = 'semi-implicit';
%   in the datafile for the fluid subproblem.
%
%   - the geometry is treated using the Geometric Convective Explicit (GCE)
%   approach; as a result, the mesh motion problem is uncoupled from the
%   fluid and solid equations. See, e.g., the following paper
%    "Crosetto, et al., Parallel algorithms for fluid-structure interaction 
%     problems in hemodynamics, SISC 2011."
%   - a condensed formulation is employed, i.e. only internal and interface
%   fluid velocity, fluid pressure, and internal solid displacement are
%   considered as degrees of freedom. See, e.g., the following paper
%    "Gee et al., Truly monolithic algebraic multigrid for fluid-structure 
%     interaction, IJNME 2011."
%   - the coupling is treated in a monolithic way
%   - the fluid equations are approximated in time by means of a
%   semi-implicit BDF scheme. Space discretization can be either P2,
%   P2Bubble (B1) or P1 with SUPG stabilization.
%   - the structure can be either linear or nonlinear; time discretization
%   is performed via the Newmark scheme. In case of nonlinearities, Newton
%   method is employed
%   ======================================================================
%
%   ======================================================================
%   ** OPT2: monolithic fully-implicit **
%
%   This options is enabled by setting 
%               data.time.nonlinearity  = 'implicit';
%   in the datafile for the fluid subproblem.
%
%   - the geometry is treated implicitly, however shape derivatives in the 
%   Jacobian matrix are not taken into account; as a result, the mesh motion
%   problem is only one-way coupled with the fluid and solid equations. 
%   - a condensed formulation is employed, i.e. only internal and interface
%   fluid velocity, fluid pressure, and internal solid displacement are
%   considered as degrees of freedom. See, e.g., the following paper
%    "Gee et al., Truly monolithic algebraic multigrid for fluid-structure 
%     interaction, IJNME 2011."
%   - the coupling is treated in a monolithic way
%   - the fluid equations are approximated in time by means of a
%   fully implicit BDF scheme. Space discretization can be either P2,
%   P2Bubble (B1) or P1 with SUPG stabilization.
%   - the structure can be either linear or nonlinear; time discretization
%   is performed via the Newmark scheme
%   - Newton method is employed to handle the nonlinearities
%   ======================================================================

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 8
    param = [];
end

if nargin < 9
    vtk_filename = [];
end


%% Read Fluid problem parameters and BCs from data_file
DATA.Fluid   = CFD_read_DataFile(data_file_F, dim, param);
if nargin < 8
    DATA.Fluid.param = [];
else
    DATA.Fluid.param = param;
end

use_SUPG = false;
if isfield(DATA.Fluid, 'Stabilization') 
    if strcmp( DATA.Fluid.Stabilization, 'SUPG' ) && strcmp(fem_F{1}, 'P1')
        use_SUPG = true;
    end
end

%% Read Solid problem parameters and BCs from data_file
DATA.Solid   = CSM_read_DataFile(data_file_S, dim, param);
if nargin < 8
    DATA.Solid.param = [];
else
    DATA.Solid.param = param;
end

%% Set quadrature order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

%% Generate Fluid and Solid mesh data structures
[ MESH.Fluid ] = buildMESH( dim, meshFluid.elements, meshFluid.vertices, ...
    meshFluid.boundaries, fem_F{1}, quad_order, DATA.Fluid, 'CFD', meshFluid.rings );
[ MESH.Solid ] = buildMESH( dim, meshSolid.elements, meshSolid.vertices, ...
    meshSolid.boundaries, fem_S, quad_order, DATA.Solid, 'CSM', meshSolid.rings );

MESH.dim  = dim;

%% Create Finite Element Spaces for the fluid velocity and pressure, solid displacement and fluid mesh displacement
[ FE_SPACE_v ] = buildFESpace( MESH.Fluid, fem_F{1}, dim, quad_order );% fluid velocity
[ FE_SPACE_p ] = buildFESpace( MESH.Fluid, fem_F{2}, 1,   quad_order );% fluid pressure
[ FE_SPACE_s ] = buildFESpace( MESH.Solid, fem_S,    dim, quad_order );% solid displacement
[ FE_SPACE_g ] = buildFESpace( MESH.Fluid, fem_F{1}, dim, quad_order );% geometry displacement

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

%% Find interface indices in the internal numbering
tmp = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
for i = 1 : dim
    tmp([FE_SPACE_v.numDofScalar*(i-1)+MESH.Fluid.dof_interface{i}]) = 1;
end

% fluid interface DoFs wrt internal numbering
MESH.Fluid.Gamma     = find(tmp(MESH.Fluid.internal_dof));
% fluid non-interface DoFs wrt internal numbering
MESH.Fluid.II        = setdiff(1:length(MESH.Fluid.internal_dof),MESH.Fluid.Gamma);
% fluid non-interface DoFs wrt original numbering
MESH.Fluid.II_global = setdiff(MESH.Fluid.internal_dof,MESH.Fluid.internal_dof(MESH.Fluid.Gamma));

tmp = zeros(FE_SPACE_s.numDof,1);
for i = 1 : dim
    tmp([FE_SPACE_s.numDofScalar*(i-1)+MESH.Solid.dof_interface{i}]) = 1;
end
% solid interface DoFs wrt internal numbering
MESH.Solid.Gamma     = find(tmp(MESH.Solid.internal_dof));
% solid non-interface DoFs wrt internal numbering
MESH.Solid.II        = setdiff(1:length(MESH.Solid.internal_dof),MESH.Solid.Gamma);
% solid non-interface DoFs wrt original numbering
MESH.Solid.II_global = setdiff(MESH.Solid.internal_dof,MESH.Solid.internal_dof(MESH.Solid.Gamma));

MESH.internal_dof    = [MESH.Fluid.II MESH.Fluid.Gamma' ...
            length(MESH.Fluid.internal_dof)+1:(length(MESH.Fluid.internal_dof)+length(MESH.Solid.II))];

% undeformed fluid mesh nodes coordinates        
Fluid_ReferenceNodes = MESH.Fluid.nodes(1:dim,:);

%% Time Setting
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
k_t  = 0;

%% Initalize Fluid Time Advance
TimeAdvanceF = BDF_TimeAdvance( BDF_orderF );

% read initial condition
v0  = [];
for k = 1 : FE_SPACE_v.numComponents
    switch dim
        case 2
            v0  = [v0; DATA.Fluid.u0{k}(  MESH.Fluid.nodes(1,:), MESH.Fluid.nodes(2,:), t0, param )'];
            
        case 3
            v0  = [v0; DATA.Fluid.u0{k}(  MESH.Fluid.nodes(1,:), MESH.Fluid.nodes(2,:), MESH.Fluid.nodes(3,:), t0, param )'];
    end
end

if isfield(DATA.Fluid, 'p0')
    switch dim
        case 2
            p0  = DATA.Fluid.p0(  MESH.Fluid.nodes(1,:), MESH.Fluid.nodes(2,:), t0, param )';
            
        case 3
            p0  = DATA.Fluid.p0(  MESH.Fluid.nodes(1,:), MESH.Fluid.nodes(2,:), MESH.Fluid.nodes(3,:), t0, param )';
    end
else
    p0 = zeros(FE_SPACE_p.numDof,1);
end

u = [v0; p0];
X_n(1:length(MESH.Fluid.internal_dof)) = u(MESH.Fluid.internal_dof);

% export initial condition (if it's the case)
if ~isempty(vtk_filename)
    CFD_export_solution(MESH.dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), ...
        MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], 0);
end

TimeAdvanceF.Initialize( v0 );
for bd = 2 : BDF_orderF
    TimeAdvanceF.Append( v0 );
end

%% Initalize Solid Time Advance
TimeAdvanceS = Newmark_TimeAdvance( DATA.Solid.time.beta, DATA.Solid.time.gamma, dt );

% read displacement and velocity initial condition
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

% export initial condition (if it's the case)
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u0, MESH.Solid.vertices, MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename,'Solid'], 0);
end

TimeAdvanceS.Initialize( u0, du0, d2u0 );
Coef_MassS = TimeAdvanceS.MassCoefficient( );


%% Initalize Geometry Time Advance
% time discretization of the ALE follows the fluid one
ALE_velocity = zeros(MESH.Fluid.numNodes*dim, 1);
d_Fn         = zeros(MESH.Fluid.numNodes*dim, 1);

% only used in the fully implicit case
TimeAdvanceG = BDF_TimeAdvance( BDF_orderF );
TimeAdvanceG.Initialize( d_Fn );
for bd = 2 : BDF_orderF
    TimeAdvanceG.Append( d_Fn );
end


%% Generate Domain Decomposition (if required)
PreconFactory = PreconditionerFactory( );
Precon        = PreconFactory.CreatePrecon(DATA.Fluid.Preconditioner.type, DATA.Fluid);

if isfield(DATA.Fluid.Preconditioner, 'type') && strcmp( DATA.Fluid.Preconditioner.type, 'AdditiveSchwarz')
    
    if isfield(DATA.Fluid.Preconditioner, 'coarse_level')
        if ~strcmp( DATA.Fluid.Preconditioner.coarse_level, 'None')
            R = FSI_overlapping_DD(MESH, DATA.Fluid.Preconditioner.num_subdomains,  ...
                     DATA.Fluid.Preconditioner.overlap_level, DATA.Fluid.Preconditioner.coarse_num_aggregates);
        else
            R = FSI_overlapping_DD(MESH, DATA.Fluid.Preconditioner.num_subdomains,  DATA.Fluid.Preconditioner.overlap_level);
        end
    else
        R = FSI_overlapping_DD(MESH, DATA.Fluid.Preconditioner.num_subdomains,  DATA.Fluid.Preconditioner.overlap_level);
    end
    
    Precon.SetRestrictions( R );
    clear R;
end

%% Assemble Solid-Extension matrix
fprintf('\n Assemble Solid Extension matrix and store LU factors ... '); 
time2          = tic; 

% create mesh motion DATA structure: young and poisson coefficients for the
% solid extension have to be set in the solid datafile
DATA.Geometry                  = DATA.Solid;
DATA.Geometry.Material_Model   = 'SEMMT';
DATA.Geometry.Stiffening_power = 0.8;

% assemble matrix
MeshMotionAssembler = CSM_Assembler( MESH.Fluid, DATA.Geometry, FE_SPACE_g );
Solid_Extension.matrix  = MeshMotionAssembler.compute_jacobian( zeros(FE_SPACE_g.numDof, 1) );

% solid-extension internal DoFs
internal_dofs_HE = [];
for k = 1 : dim
    internal_dofs_HE       = [ internal_dofs_HE (k-1)*FE_SPACE_v.numDofScalar+setdiff(1:FE_SPACE_v.numDofScalar, ...
        [MESH.Fluid.dof_interface{k}; MESH.ALE_dirichlet{k}])];
end
Solid_Extension.internal_dofs = internal_dofs_HE;

% compute LU factorization and store it
[Solid_Extension.L , Solid_Extension.U,...
    Solid_Extension.perm , q ]   = lu(Solid_Extension.matrix(internal_dofs_HE,internal_dofs_HE), 'vector');
Solid_Extension.invp             = 0*q ;
Solid_Extension.invp(q)          = 1:length(q);

time2          = toc(time2);
fprintf('%f s \n',time2);
fprintf('-------------------------------------------\n');

%% Interface transfer matrices: solid to fluid and viceversa

% empty matrices for zero matrix blocks in the monolithic system
Z_FS  = sparse(length(MESH.Fluid.II), length(MESH.Solid.II));
Z_SF  = sparse(length(MESH.Solid.II), length(MESH.Fluid.II));

% transfer matrix from solid to fluid
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

%% Create Solid Assembler Object
SolidModel = CSM_Assembler( MESH.Solid, DATA.Solid, FE_SPACE_s );

% Assemble mass matrix
fprintf('\n Assembling S-mass matrix... ');
t_assembly = tic;
M_s    =  SolidModel.compute_mass();
M_s    =  M_s * DATA.Solid.Density;
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

% if the material model is linear elasticity, assemble stiffness matrix
if strcmp( DATA.Solid.Material_Model, 'Linear' )
    A_s = SolidModel.compute_jacobian( zeros(FE_SPACE_s.numDof, 1) );
end

% Assemble Robin BC (if it's the case)
A_robin = SolidModel.assemble_ElasticRobinBC();

%% Initialize Linear Solver
LinSolver  = LinearSolver( DATA.Fluid.LinearSolver );

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

%% PreProcessing for Boundary Flow Rates computations
compute_FlowRates = 0 ;
if isfield(DATA.Fluid, 'Output') && isfield(DATA.Fluid.Output, 'FlowRates')
    if DATA.Fluid.Output.FlowRates.computeFlowRates == 1
        compute_FlowRates = true;
    end
end

if compute_FlowRates
    
    fileFlowRates = fopen(DATA.Fluid.Output.FlowRates.filename, 'w+');
    fprintf(fileFlowRates, 'Time');
    
    for l = 1 : length(DATA.Fluid.Output.FlowRates.flag)
        fprintf(fileFlowRates, '         Flag %d', DATA.Fluid.Output.FlowRates.flag(l) );
    end
    fprintf(fileFlowRates, '         Sum');
    
    fprintf(fileFlowRates, '\n%1.3e', t);
    for l = 1 : length(DATA.Fluid.Output.FlowRates.flag)
        FlowRate(l)  = CFD_computeFlowRate(u, MESH.Fluid, FE_SPACE_v, FE_SPACE_p, DATA.Fluid.Output.FlowRates.flag(l));
        fprintf(fileFlowRates, '    %1.3e', FlowRate(l) );
    end
    fprintf(fileFlowRates, '    %1.3e', sum( FlowRate ) );
    
end

%% Load Prestress data
use_Prestress = false ;
if isfield(DATA.Solid, 'Prestress')
    if DATA.Solid.Prestress == true
        use_Prestress = true;
    end
end

if use_Prestress

	load R_P R_P;
	load J_P J_P;
    fprintf('\n\n -- Prestress loaded -- \n');

else
	R_P = zeros(FE_SPACE_s.numDof,1);
	J_P = sparse(FE_SPACE_s.numDof,FE_SPACE_s.numDof);
end

X_nk = X_n;

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
    alphaF = TimeAdvanceF.GetCoefficientDerivative();
    
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
            C_NS = alphaF/dt * M + A_Stokes + C1;
            
            % Assemble SUPG contributes
            if use_SUPG
                fprintf('\n   -- Fluid_Assembling SUPG Terms ... ');
                t_assembly = tic;
                [A_SUPG, F_SUPG] = FluidModel.compute_SUPG_semiimplicit( v_extrapolated - ALE_velocity, v_BDF, dt, alphaF);
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
                    d_alpha = TimeAdvanceS.M_dU - dt * DATA.Solid.time.gamma * Csi + dt* (1-DATA.Solid.time.gamma)*TimeAdvanceS.M_d2U;

                    F_L     = d_alpha(MESH.Solid.internal_dof(MESH.Solid.Gamma));
                    alpha   = DATA.Solid.time.gamma / (dt * DATA.Solid.time.beta);
                    
                    % Interface Solid Stiffness expressed in the fluid numbering
                    S_GG     = IdGamma_SF * (C_STR(MESH.Solid.Gamma,MESH.Solid.Gamma) * IdGamma_FS);

                    % Interface/Internal Solid Stiffness expressed in fluid/solid numbering
                    S_GI     = IdGamma_SF *  C_STR(MESH.Solid.Gamma,MESH.Solid.II);
                    
                    % Form Monolothic System
                    FSI_M    = [C_NS_in(MESH.Fluid.II,MESH.Fluid.II)                      C_NS_in(MESH.Fluid.II,MESH.Fluid.Gamma)            Z_FS ;...
                        C_NS_in(MESH.Fluid.Gamma,MESH.Fluid.II)      C_NS_in(MESH.Fluid.Gamma,MESH.Fluid.Gamma)+1/alpha*S_GG                 S_GI  ;...
                        Z_SF                                                    1/alpha*C_STR(MESH.Solid.II,MESH.Solid.Gamma)*IdGamma_FS     C_STR(MESH.Solid.II,MESH.Solid.II)];
                    
                    FSI_F    = [F_NS_in(MESH.Fluid.II); ...
                        F_NS_in(MESH.Fluid.Gamma)+IdGamma_SF*F_S(MESH.Solid.Gamma)+1/alpha*(IdGamma_SF*(C_STR(MESH.Solid.Gamma,MESH.Solid.Gamma)*F_L)); ...
                        F_S(MESH.Solid.II)+1/alpha*C_STR(MESH.Solid.II,MESH.Solid.Gamma)*F_L];
                                        
                    % Solve Monolothic System
                    fprintf('\n -- Solve A x = b ... ');
                    Precon.Build( FSI_M );
                    fprintf('\n      time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
                    LinSolver.SetPreconditioner( Precon );
                    X_nk(MESH.internal_dof) = LinSolver.Solve( FSI_M, FSI_F, X_n(MESH.internal_dof) );
                    fprintf('\n      time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
                    
                otherwise
                    
                    [~, ~, DisplacementDir_np1] = CSM_ApplyBC([], [], FE_SPACE_s, MESH.Solid, DATA.Solid, t);
                    
                    Csi           = TimeAdvanceS.RhsContribute( );
                    F_ext         = SolidModel.compute_volumetric_forces( t );
                    F_S           = F_ext + M_s * Csi;
                    
                    % Get displacment d^alpha on the interface
                    d_alpha = TimeAdvanceS.M_dU - dt * DATA.Solid.time.gamma * Csi ...
                                   + dt* (1-DATA.Solid.time.gamma)*TimeAdvanceS.M_d2U;
                    F_L     = d_alpha(MESH.Solid.internal_dof(MESH.Solid.Gamma));
                    
                    alpha   = DATA.Solid.time.gamma / (dt * DATA.Solid.time.beta);
                    
                    % Prepare for Newton's method
                    X_nk                          = X_n;
                    dX                            = 0*X_n;
                    
                    % fluid current iteration
                    u_nk                           = 0*u;
                    u_nk(MESH.Fluid.internal_dof)  = X_nk(1:length(MESH.Fluid.internal_dof));
                    u_nk(MESH.Fluid.Dirichlet_dof) = v_D;
                    
                    % solid current iteration
                    d_nk                               = TimeAdvanceS.M_U;
                    d_nk(MESH.Solid.Dirichlet_dof)     = DisplacementDir_np1;
                    
                    k                         = 1;
                    norm_k                    = tol + 1;
                    
                    % Start Newton's method
                    while( k <= maxIter && (norm_k > tol || isnan(norm_k) ) )
                        
                        % Assemble Solid Tangent matrix and internal forces vector
                        fprintf('\n   -- Solid_Assembling Residual and Jacobian ... ');
                        t_assembly = tic;
                        dA = SolidModel.compute_jacobian(  d_nk  );
                        GS = SolidModel.compute_internal_forces( d_nk );
                        t_assembly = toc(t_assembly);
                        fprintf('done in %3.3f s\n', t_assembly);
                        
                        dG_STR    = Coef_MassS * M_s + dA;
                        G_S       = Coef_MassS * M_s * d_nk + GS - F_S;
                        
                        % Apply Solid boundary conditions
                        [dG_STR, G_S] = CSM_ApplyBC(dG_STR, -G_S, FE_SPACE_s, MESH.Solid, DATA.Solid, t, 1);
                        
                        % Apply Fluid boundary conditions
                        [dG_NS, G_NS]   =  CFD_ApplyBC(C_NS, -(C_NS*u_nk - F_NS), FE_SPACE_v, FE_SPACE_p, MESH.Fluid, DATA.Fluid, t, 1);
                        
                        % Interface Solid Stiffness expressed in the fluid numbering
                        S_GG     = IdGamma_SF * (dG_STR(MESH.Solid.Gamma,MESH.Solid.Gamma) * IdGamma_FS);
                        
                        % Interface/Internal Solid Stiffness expressed in fluid/solid numbering
                        S_GI     = IdGamma_SF *  dG_STR(MESH.Solid.Gamma,MESH.Solid.II);
                        F_L2     = F_L - IdGamma_FS*X_nk(MESH.Fluid.Gamma) + alpha*d_nk(MESH.Solid.internal_dof(MESH.Solid.Gamma));
                        
                        % Form Monolothic System
                        dG_FSI    = [dG_NS(MESH.Fluid.II,MESH.Fluid.II)                      dG_NS(MESH.Fluid.II,MESH.Fluid.Gamma)                    Z_FS ;...
                            dG_NS(MESH.Fluid.Gamma,MESH.Fluid.II)      dG_NS(MESH.Fluid.Gamma,MESH.Fluid.Gamma)+1/alpha*S_GG               S_GI  ;...
                            Z_SF                                                    1/alpha*dG_STR(MESH.Solid.II,MESH.Solid.Gamma)*IdGamma_FS   dG_STR(MESH.Solid.II,MESH.Solid.II)   ];
                        
                        G_FSI    = [G_NS(MESH.Fluid.II); ...
                            G_NS(MESH.Fluid.Gamma)+IdGamma_SF*G_S(MESH.Solid.Gamma)+1/alpha*(IdGamma_SF*(dG_STR(MESH.Solid.Gamma,MESH.Solid.Gamma)*F_L2)); ...
                            G_S(MESH.Solid.II)+1/alpha*dG_STR(MESH.Solid.II,MESH.Solid.Gamma)*F_L2];
                        
                        % Solve Monolothic System
                        fprintf('\n -- Solve J x = -R ... ');
                        Precon.Build( dG_FSI );
                        fprintf('\n      time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
                        LinSolver.SetPreconditioner( Precon );
                        dX(MESH.internal_dof) = LinSolver.Solve( dG_FSI, G_FSI );
                        fprintf('\n      time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
                        
                        % Update solution
                        norm_k   = norm(dX)/norm(X_nk);
                        X_nk     = X_nk + dX;
                        
                        u_nk(MESH.Fluid.internal_dof)             = X_nk(1:length(MESH.Fluid.internal_dof));
                        u_nk(MESH.Fluid.Dirichlet_dof) = v_D;
                        
                        d_nk(MESH.Solid.II_global)              = X_nk(1+length(MESH.Fluid.internal_dof):end);
                        d_nk(MESH.Solid.Dirichlet_dof)          = DisplacementDir_np1;
                        d_nk(MESH.Solid.internal_dof(MESH.Solid.Gamma))     = 1/alpha * ( IdGamma_FS*X_nk(MESH.Fluid.Gamma) -  F_L);
                        
                        fprintf('\n   Iteration  k= %d;  norm(dX)/norm(Xk) = %1.2e, normRES = %1.2e \n',k,full(norm_k), full(norm(G_FSI)))
                        k = k + 1;
                        
                    end
            end
            
        case 'implicit'
            
            % get lifting for F-velocity and S_displacement
            [~, ~, v_D]                 = CFD_ApplyBC([], [], FE_SPACE_v, FE_SPACE_p, MESH.Fluid, DATA.Fluid, t);            
            [~, ~, DisplacementDir_np1] = CSM_ApplyBC([], [], FE_SPACE_s, MESH.Solid, DATA.Solid, t);
            
            Csi           = TimeAdvanceS.RhsContribute( );
            F_ext         = SolidModel.compute_volumetric_forces( t );
            F_S           = F_ext + M_s * Csi;
            
            % Get displacment d^alpha on the interface
            d_alpha = TimeAdvanceS.M_dU - dt * DATA.Solid.time.gamma * Csi ...
                + dt* (1-DATA.Solid.time.gamma)*TimeAdvanceS.M_d2U;
            F_L     = d_alpha(MESH.Solid.internal_dof(MESH.Solid.Gamma));
            
            alpha   = DATA.Solid.time.gamma / (dt * DATA.Solid.time.beta);
            
            % Prepare for Newton's method
            X_nk                          = X_n;
            dX                            = 0*X_n;
            
            % fluid current iteration
            u_nk                           = 0*u;
            u_nk(MESH.Fluid.internal_dof)  = X_nk(1:length(MESH.Fluid.internal_dof));
            u_nk(MESH.Fluid.Dirichlet_dof) = v_D;
            
            % solid current iteration
            d_nk                               = TimeAdvanceS.M_U;
            d_nk(MESH.Solid.Dirichlet_dof)     = DisplacementDir_np1;
            
            ALE_rhsBDF = TimeAdvanceG.RhsContribute( );
            
            k        = 1;
            norm_k   = tol + 1;
            
            % Start Newton's method
            while( k <= maxIter && (norm_k > tol || isnan(norm_k) ) )
                
                % Assemble Solid Tangent matrix and internal forces vector
                fprintf('\n   -- Solid_Assembling Residual and Jacobian ... ');
                t_assembly = tic;
                dA = SolidModel.compute_jacobian(  d_nk  );
                GS = SolidModel.compute_internal_forces( d_nk );
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                
                dG_STR    = Coef_MassS * M_s + dA + A_robin + J_P;
                G_S       = Coef_MassS * M_s * d_nk + GS + A_robin * d_nk + R_P + J_P * d_nk - F_S;

                % Apply Solid boundary conditions
                [dG_STR, G_S] = CSM_ApplyBC(dG_STR, -G_S, FE_SPACE_s, MESH.Solid, DATA.Solid, t, 1);
                
                % Update Fluid Matrices
                FluidModel = CFD_Assembler( MESH.Fluid, DATA.Fluid, FE_SPACE_v, FE_SPACE_p );
                
                fprintf('\n   -- Fluid_Assembling volumetric forces... ');
                t_assembly = tic;
                F_gravity = FluidModel.compute_external_forces(t);
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                
                fprintf('\n   -- Fluid_Assembling Stokes terms... ');
                t_assembly = tic;
                [A_Stokes] = FluidModel.compute_Stokes_matrix();
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                
                fprintf('\n   -- Fluid_Assembling Mass matrix... ');
                t_assembly = tic;
                Mv = FluidModel.compute_mass_velocity();
                Mp = FluidModel.compute_mass_pressure();
                M  = blkdiag(DATA.Fluid.density * Mv, 0*Mp);
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                                
                fprintf('\n   -- Fluid_Assembling Convective Term... ');
                t_assembly = tic;
                [C1, C2]    = FluidModel.compute_convective_matrix_ALE( u_nk(1:FE_SPACE_v.numDof),  ALE_velocity);
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                
                F_NS = 1/dt * M * (alphaF*u_nk - u_BDF) + A_Stokes * u_nk + C1 * u_nk - F_gravity;
                C_NS = alphaF/dt * M + A_Stokes + C1 + C2;
                
                % Assemble SUPG contributes
                if use_SUPG
                    fprintf('\n   -- Fluid_Assembling SUPG Terms ... ');
                    t_assembly = tic;
                    [A_SUPG, F_SUPG] = FluidModel.compute_SUPG_implicit_ALE( u_nk, ALE_velocity, v_BDF, dt, alphaF);
                    t_assembly = toc(t_assembly);
                    fprintf('done in %3.3f s\n', t_assembly);
                    
                    C_NS             = C_NS + A_SUPG;
                    F_NS             = F_NS + F_SUPG;
                end
                
                % Apply Fluid boundary conditions
                [dG_NS, G_NS]   =  CFD_ApplyBC(C_NS, -F_NS, FE_SPACE_v, FE_SPACE_p, MESH.Fluid, DATA.Fluid, t, 1, u);
                
                fprintf('\n   -- Form monolithic system ... ');
                t_assembly = tic;
                % Interface Solid Stiffness expressed in the fluid numbering
                S_GG     = IdGamma_SF * (dG_STR(MESH.Solid.Gamma,MESH.Solid.Gamma) * IdGamma_FS);
                
                % Interface/Internal Solid Stiffness expressed in fluid/solid numbering
                S_GI     = IdGamma_SF *  dG_STR(MESH.Solid.Gamma,MESH.Solid.II);
                F_L2     = F_L - IdGamma_FS*X_nk(MESH.Fluid.Gamma) + alpha*d_nk(MESH.Solid.internal_dof(MESH.Solid.Gamma));
                
                % Form Monolothic System
                dG_FSI    = [dG_NS(MESH.Fluid.II,MESH.Fluid.II)                      dG_NS(MESH.Fluid.II,MESH.Fluid.Gamma)                    Z_FS ;...
                    dG_NS(MESH.Fluid.Gamma,MESH.Fluid.II)      dG_NS(MESH.Fluid.Gamma,MESH.Fluid.Gamma)+1/alpha*S_GG               S_GI  ;...
                    Z_SF                                                    1/alpha*dG_STR(MESH.Solid.II,MESH.Solid.Gamma)*IdGamma_FS   dG_STR(MESH.Solid.II,MESH.Solid.II)   ];
                
                G_FSI    = [G_NS(MESH.Fluid.II); ...
                    G_NS(MESH.Fluid.Gamma)+IdGamma_SF*G_S(MESH.Solid.Gamma)+1/alpha*(IdGamma_SF*(dG_STR(MESH.Solid.Gamma,MESH.Solid.Gamma)*F_L2)); ...
                    G_S(MESH.Solid.II)+1/alpha*dG_STR(MESH.Solid.II,MESH.Solid.Gamma)*F_L2];
                
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                    
                % Solve Monolothic System
                fprintf('\n -- Solve J x = -R ... ');
                Precon.Build( dG_FSI );
                fprintf('\n      time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
                LinSolver.SetPreconditioner( Precon );
                dX(MESH.internal_dof) = LinSolver.Solve( dG_FSI, G_FSI );
                fprintf('\n      time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
                
                % Update solution
                norm_k   = norm(dX)/norm(X_nk);
                X_nk     = X_nk + dX;
                
                % Update current velocity and pressure
                u_nk(MESH.Fluid.internal_dof)             = X_nk(1:length(MESH.Fluid.internal_dof));
                u_nk(MESH.Fluid.Dirichlet_dof) = v_D;
                
                % Update current solid displacement
                d_nk(MESH.Solid.II_global)              = X_nk(1+length(MESH.Fluid.internal_dof):end);
                d_nk(MESH.Solid.Dirichlet_dof)          = DisplacementDir_np1;
                d_nk(MESH.Solid.internal_dof(MESH.Solid.Gamma))     = 1/alpha * ( IdGamma_FS*X_nk(MESH.Fluid.Gamma) -  F_L);
                
                % Deform Fluid mesh by Solid-Extension Mesh Motion technique
                d_F = FSI_SolidExtension(MESH, d_nk, Solid_Extension);
                Fluid_def_nodes    = Fluid_ReferenceNodes + d_F;
                Fluid_def_vertices = Fluid_def_nodes(1:dim, 1:MESH.Fluid.numVertices);
                
                d_F      = reshape(d_F',dim*MESH.Fluid.numNodes,1);
                
                % Compute Fluid mesh velocity by BDF formula
                ALE_velocity  =  1/dt * ( alphaF*d_F - ALE_rhsBDF );
                
                % Update Fluid MESH jacobian etc
                MESH.Fluid.vertices = Fluid_def_vertices;
                MESH.Fluid.nodes    = Fluid_def_nodes;
                [MESH.Fluid.jac, MESH.Fluid.invjac, MESH.Fluid.h] = geotrasf(dim, MESH.Fluid.vertices, MESH.Fluid.elements);
                
                fprintf('\n   Iteration  k= %d;  norm(dX)/norm(Xk) = %1.2e, normRES = %1.2e \n',k,full(norm_k), full(norm(G_FSI)))
                k = k + 1;
                
            end
                
    end
    
    norm_n   = norm(X_nk - X_n) / norm(X_n);
    X_n      = X_nk;

    %% Export solid displacement on reference mesh
    Displacement_np1                                   = zeros(FE_SPACE_s.numDof,1);
    Displacement_np1(MESH.Solid.II_global)             = X_n(1+length(MESH.Fluid.internal_dof):end);
    Displacement_np1(MESH.Solid.Dirichlet_dof)         = DisplacementDir_np1;
    Displacement_np1(MESH.Solid.internal_dof(MESH.Solid.Gamma))   = 1/alpha * ( IdGamma_FS*X_n(MESH.Fluid.Gamma) -  F_L);
    
    % Export to VTK
    if ~isempty(vtk_filename)
        CSM_export_solution(MESH.dim, Displacement_np1, MESH.Solid.vertices, ...
            MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename, 'Solid'], k_t);
    end
    
    % update time advance
    TimeAdvanceS.Update( Displacement_np1 );
    
    %% If GCE is used, update Fluid Mesh
    if strcmp(DATA.Fluid.time.nonlinearity, 'semi-implicit')
        
        % Deform Fluid mesh by Solid-Extension Mesh Motion technique
        d_F = FSI_SolidExtension(MESH, Displacement_np1, Solid_Extension);
        Fluid_def_nodes    = Fluid_ReferenceNodes + d_F;
        Fluid_def_vertices = Fluid_def_nodes(1:dim, 1:MESH.Fluid.numVertices);
        
        d_F      = reshape(d_F',dim*MESH.Fluid.numNodes,1);
        
        % Compute Fluid mesh velocity: w = 1/dt * ( d_f^(n+1) - d_f^n )
        % This part should be checked! inconsistent with BDF integrator
        ALE_velocity  =  1/dt * ( d_F - d_Fn );
        d_Fn          =  d_F;
        
        % Update Fluid MESH
        MESH.Fluid.vertices = Fluid_def_vertices;
        MESH.Fluid.nodes    = Fluid_def_nodes;
        % update mesh jacobian, determinant and inverse
        [MESH.Fluid.jac, MESH.Fluid.invjac, MESH.Fluid.h] = geotrasf(dim, MESH.Fluid.vertices, MESH.Fluid.elements);
    end
    
    % Update geometry time advance
    TimeAdvanceG.Append( d_F );
    d_Fn          =  d_F;

    %% Export Fluid velocity and pressure on deformed mesh
    u(MESH.Fluid.internal_dof)  = X_n(1:length(MESH.Fluid.internal_dof));
    u(MESH.Fluid.Dirichlet_dof) = v_D;
        
    % Export to VTK
    if ~isempty(vtk_filename)
        CFD_export_solution(dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), ...
            MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], k_t);
    end
    
    % Update fluid time advance
    TimeAdvanceF.Append( u(1:FE_SPACE_v.numDof) );
    
    %% Compute Aerodynamic Forces
    if compute_AerodynamicForces
        
        if strcmp(DATA.Fluid.time.nonlinearity,'implicit')
            C_NS = 0*C_NS;
            F_NS = -F_NS;
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

    if compute_FlowRates        
        fprintf(fileFlowRates, '\n%1.3e', t);
        for l = 1 : length(DATA.Fluid.Output.FlowRates.flag)
            FlowRate(l)  = CFD_computeFlowRate(u, MESH.Fluid, FE_SPACE_v, FE_SPACE_p, DATA.Fluid.Output.FlowRates.flag(l));
            fprintf(fileFlowRates, '    %1.3e', FlowRate(l) );
        end
        fprintf(fileFlowRates, '    %1.3e', sum( FlowRate ) );
    end

    iter_time = toc(iter_time);
    fprintf('\nnorm(U^(n+1) - U^(n))/norm(U^(n)) = %2.3e -- Iteration time: %3.2f s \n',full(norm_n),iter_time);
    
    X = [u; Displacement_np1];
    
end

fprintf('\n************************************************************************* \n');

return
