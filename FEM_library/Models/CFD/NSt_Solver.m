function [u, MESH, DATA] = NSt_Solver(dim, elements, vertices, boundaries, fem, data_file, param, vtk_filename)
%NST_SOLVER unsteady Navier-Stokes Equations solver

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

%% Gather Time Setting
BDF_order = DATA.time.BDF_order;
t0        = DATA.time.t0;
dt        = DATA.time.dt;
tf        = DATA.time.tf;
t         = DATA.time.t0;
k_t       = 0;

BDFhandler = BDF_TimeAdvance( BDF_order);

v0  = [];
for k = 1 : FE_SPACE_v.numComponents
    switch dim
        case 2
            v0  = [v0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), t0, param )'];
            
        case 3
            v0  = [v0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];
    end
end
u = [v0; zeros(FE_SPACE_p.numDof,1)];
if ~isempty(vtk_filename)
    CFD_export_solution(MESH.dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, 0);
end

BDFhandler.Initialize( v0 );
for bd = 2 : BDF_order
    BDFhandler.Append( v0 );
end

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',length(MESH.internal_dof));
fprintf(' * Number of timesteps =  %d\n', (tf-t0)/dt);
fprintf(' * BDF Order           =  %d\n', BDF_order);
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

fprintf('\n Assembling mass matrix... ');
t_assembly = tic;
Mv = FluidModel.compute_mass_velocity();
Mp = FluidModel.compute_mass_pressure();
M  = blkdiag(DATA.density * Mv, 0*Mp);
t_assembly = toc(t_assembly);
fprintf('done in %3.3f s', t_assembly);

%% Initialize Linear Solver
LinSolver = LinearSolver( DATA.LinearSolver );


%% PreProcessing for Drag and Lift Computation
compute_AerodynamicForces = 0 ;
if isfield(DATA, 'Output') && isfield(DATA.Output, 'DragLift')
    if DATA.Output.DragLift.computeDragLift == 1
        compute_AerodynamicForces = true;
    end
end

if compute_AerodynamicForces
    AeroF_x(k_t+1)  = 0;
    AeroF_y(k_t+1)  = 0;
    AeroF_z(k_t+1)  = 0;
    dofs_drag    = [];
    
    for j = 1 : length(DATA.Output.DragLift.flag)
        Dirichlet_side         = find(MESH.boundaries(MESH.bc_flag_row,:) == DATA.Output.DragLift.flag(j));
        Dirichlet_side         = unique(Dirichlet_side);
        Dirichlet_dof          = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
        dofs_drag              = [dofs_drag; Dirichlet_dof(:)];
    end
    dofs_drag = unique(dofs_drag);
    
    fileDragLift = fopen(DATA.Output.DragLift.filename, 'w+');
    fprintf(fileDragLift, 'Time          F_x          F_y          F_z');
    fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
end

%% PreProcessing for Boundary Flow Rates computations
compute_FlowRates = 0 ;
if isfield(DATA, 'Output') && isfield(DATA.Output, 'FlowRates')
    if DATA.Output.FlowRates.computeFlowRates == 1
        compute_FlowRates = true;
    end
end

if compute_FlowRates
    
    fileFlowRates = fopen(DATA.Output.FlowRates.filename, 'w+');
    fprintf(fileFlowRates, 'Time');
    
    for l = 1 : length(DATA.Output.FlowRates.flag)
        fprintf(fileFlowRates, '         Flag %d', DATA.Output.FlowRates.flag(l) );
    end
    fprintf(fileFlowRates, '         Sum');
    
    fprintf(fileFlowRates, '\n%1.3e', t);
    for l = 1 : length(DATA.Output.FlowRates.flag)
        FlowRate(l)  = CFD_computeFlowRate(u, MESH, FE_SPACE_v, FE_SPACE_p, DATA.Output.FlowRates.flag(l));
        fprintf(fileFlowRates, '    %1.3e', FlowRate(l) );
    end
    fprintf(fileFlowRates, '    %1.3e', sum( FlowRate ) );
    
end

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n',t0,t,tf);
    
    v_BDF = BDFhandler.RhsContribute( );
    u_BDF = [v_BDF; zeros(FE_SPACE_p.numDof,1)];
    alpha = BDFhandler.GetCoefficientDerivative();
    
    switch DATA.time.nonlinearity
        
        case 'semi-implicit'
            
            v_extrapolated = BDFhandler.Extrapolate();
            U_k            = zeros(totSize,1);
            
            % Assemble matrix and right-hand side
            fprintf('\n -- Assembling Convective Term... ');
            t_assembly = tic;
            [C1] = FluidModel.compute_convective_Oseen_matrix( v_extrapolated );
            t_assembly = toc(t_assembly);
            fprintf('done in %3.3f s\n', t_assembly);
            
            F_NS = 1/dt * M * u_BDF;
            C_NS = alpha/dt * M + A_Stokes + C1;
            
            if use_SUPG
                
                fprintf('\n -- Assembling SUPG Terms ... ');
                t_assembly = tic;
                [A_SUPG, F_SUPG] = FluidModel.compute_SUPG_semiimplicit( v_extrapolated, v_BDF, dt, alpha);
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                
                C_NS             = C_NS + A_SUPG;
                F_NS             = F_NS - F_SUPG;
            end
            
            % Apply boundary conditions
            fprintf('\n -- Apply boundary conditions ... ');
            t_assembly = tic;
            [A, b, u_D]   =  CFD_ApplyBC(C_NS, F_NS, FE_SPACE_v, FE_SPACE_p, MESH, DATA, t, 0, u);
            t_assembly = toc(t_assembly);
            fprintf('done in %3.3f s\n', t_assembly);
            
            % Solve
            fprintf('\n -- Solve A x = b ... ');
            Precon.Build( A );
            fprintf('\n      time to build the preconditioner %3.3f s \n', Precon.GetBuildTime());
            LinSolver.SetPreconditioner( Precon );
            U_k(MESH.internal_dof) = LinSolver.Solve( A, b, u );
            fprintf('\n      time to solve the linear system in %3.3f s \n', LinSolver.GetSolveTime());
            U_k(MESH.Dirichlet_dof) = u_D;
            
            fprintf('\n -- Norm(U_np1 - U_n) / Norm( U_n ) = %1.2e \n', norm(U_k - u) / norm(u));
            
        case 'implicit'
            
            % Nonlinear Iterations
            tol        = DATA.NonLinearSolver.tol;
            resRelNorm = tol + 1;
            incrNorm   = tol + 1;
            maxIter    = DATA.NonLinearSolver.maxit;
            k          = 1;
            
            [~, ~, u_D]   =  CFD_ApplyBC([], [], FE_SPACE_v, FE_SPACE_p, MESH, DATA, t);
            dU             = zeros(totSize,1);
            U_k            = u;
            U_k(MESH.Dirichlet_dof) = u_D;
            
            % Assemble matrix and right-hand side
            fprintf('\n -- Assembling Convective terms... ');
            t_assembly = tic;
            [C1, C2] = FluidModel.compute_convective_matrix( U_k );
            t_assembly = toc(t_assembly);
            fprintf('done in %3.3f s\n', t_assembly);
            
            Residual = 1/dt * M * (alpha*U_k - u_BDF) + A_Stokes * U_k + C1 * U_k;
            Jacobian = alpha/dt * M + A_Stokes + C1 + C2;
            
            if use_SUPG
                fprintf('\n -- Assembling SUPG Terms ... ');
                t_assembly = tic;
                [A_SUPG, F_SUPG] = FluidModel.compute_SUPG_implicit( U_k, v_BDF, dt, alpha);
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                
                Jacobian    = Jacobian + A_SUPG;
                Residual    = Residual + F_SUPG;
            end
            
            % Apply boundary conditions
            fprintf('\n -- Apply boundary conditions ... ');
            t_assembly = tic;
            [A, b]   =  CFD_ApplyBC(Jacobian, -Residual, FE_SPACE_v, FE_SPACE_p, MESH, DATA, t, 1, u);
            t_assembly = toc(t_assembly);
            fprintf('done in %3.3f s\n', t_assembly);
            
            res0Norm = norm(b);
            
            fprintf('\n============ Start Newton Iterations ============\n\n');
            while (k <= maxIter && (incrNorm > tol || resRelNorm > tol))
                
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
                
                Residual = 1/dt * M * (alpha*U_k - u_BDF) + A_Stokes * U_k + C1 * U_k;
                Jacobian = alpha/dt * M + A_Stokes + C1 + C2;
                
                if use_SUPG
                    
                    fprintf('\n   -- Assembling SUPG Terms ... ');
                    t_assembly = tic;
                    [A_SUPG, F_SUPG] = FluidModel.compute_SUPG_implicit( U_k, v_BDF, dt, alpha);
                    t_assembly = toc(t_assembly);
                    fprintf('done in %3.3f s\n', t_assembly);
                    
                    Jacobian    = Jacobian + A_SUPG;
                    Residual    = Residual + F_SUPG;
                end
                
                % Apply boundary conditions
                fprintf('\n   -- Apply boundary conditions ... ');
                t_assembly = tic;
                [A, b]   =  CFD_ApplyBC(Jacobian, -Residual, FE_SPACE_v, FE_SPACE_p, MESH, DATA, t, 1, u);
                t_assembly = toc(t_assembly);
                fprintf('done in %3.3f s\n', t_assembly);
                
                resRelNorm = norm(b) / res0Norm;
                
                fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k, full(incrNorm), full(resRelNorm));
                k = k + 1;
                
            end
            fprintf('\n -- Norm(U_np1 - U_n) / Norm( U_n ) = %1.2e \n', norm(U_k - u) / norm(u));

    end
    u = U_k;
    
    %% Update BDF
    BDFhandler.Append( u(1:FE_SPACE_v.numDof) );
   
    %% Export to VTK
    if ~isempty(vtk_filename)
        CFD_export_solution(MESH.dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, k_t);
    end
       
    %% Compute_DragLift
    if compute_AerodynamicForces
        
        if strcmp(DATA.time.nonlinearity,'implicit')
            C_NS = 0*Jacobian;
            F_NS = -Residual;
        end
        Z              = zeros(FE_SPACE_v.numDofScalar,1);
        Z(dofs_drag)   = 1;
        
        W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
        W(1:FE_SPACE_v.numDofScalar)        = Z;
        AeroF_x(k_t+1) = DATA.Output.DragLift.factor*(W'*(-C_NS*u + F_NS));
        
        W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
        W(FE_SPACE_v.numDofScalar+[1:FE_SPACE_v.numDofScalar])  = Z;
        AeroF_y(k_t+1)  = DATA.Output.DragLift.factor*(W'*(-C_NS*u  + F_NS));
        
        if MESH.dim == 3
            W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
            W(2*FE_SPACE_v.numDofScalar+[1:FE_SPACE_v.numDofScalar])  = Z;
            AeroF_z(k_t+1)  = DATA.Output.DragLift.factor*(W'*(-C_NS*u  + F_NS));
        else
            AeroF_z(k_t+1) = 0.0;
        end
        
        fprintf('\n *** F_x = %e, F_y = %e, F_z = %e *** \n',  AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
        fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
    end
    
    if compute_FlowRates
        fprintf(fileFlowRates, '\n%1.3e', t);
        for l = 1 : length(DATA.Output.FlowRates.flag)
            FlowRate(l)  = CFD_computeFlowRate(u, MESH, FE_SPACE_v, FE_SPACE_p, DATA.Output.FlowRates.flag(l));
            fprintf(fileFlowRates, '    %1.3e', FlowRate(l) );
        end
        fprintf(fileFlowRates, '    %1.3e', sum( FlowRate ) );
    end

    iter_time = toc(iter_time);
    fprintf('\n-------------- Iteration time: %3.2f s -----------------',iter_time);
    
end

if compute_AerodynamicForces
    fclose(fileDragLift);
end

if compute_FlowRates
    fclose(fileFlowRates);
end

%% Compute and save Fluid Load on the boundary (for prestressing in FSI simulations)
export_FluidLoad = 0 ;
if isfield(DATA, 'Output') && isfield(DATA.Output, 'ExportFinalLoad')
    if DATA.Output.ExportFinalLoad.computeLoad == 1
        export_FluidLoad = true;
    end
end

if export_FluidLoad
    dofs_load    = [];
    
    for j = 1 : length(DATA.Output.ExportFinalLoad.flag)
        load_side         = find(MESH.boundaries(MESH.bc_flag_row,:) == DATA.Output.ExportFinalLoad.flag(j));
        load_side         = unique(load_side);
        load_dof          = MESH.boundaries(1:MESH.numBoundaryDof,load_side);
        dofs_load         = [dofs_load; load_dof(:)];
    end
    dofs_load = unique(dofs_load);
    
    if strcmp(DATA.time.nonlinearity,'implicit')
        F_NS = Residual;
    end
    FluidLoad = F_NS;%(dofs_load);
    save(DATA.Output.ExportFinalLoad.filename, 'FluidLoad');
end

fprintf('\n************************************************************************* \n');

return
