%% PHYSICAL DATA

% Inlet inward normal
N1  = [0.3229 0.8899 0.3209];

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y);
data.force{2} = @(x, y, z, t, param)(0.*x.*y);
data.force{3} = @(x, y, z, t, param)(0.*x.*y);

% Dirichlet
data.bcDir{1,5} = @(x, y, z, t, param)( Flow_Rate( t )  * N1(1) * 1 + 0.*x.*y); 
data.bcDir{2,5} = @(x, y, z, t, param)( Flow_Rate( t )  * N1(2) * 1 + 0.*x.*y); 
data.bcDir{3,5} = @(x, y, z, t, param)( Flow_Rate( t )  * N1(3) * 1 + 0.*x.*y); 

data.bcDir{1,200} = @(x, y, z, t, param)(0.*x.*y); 
data.bcDir{2,200} = @(x, y, z, t, param)(0.*x.*y); 
data.bcDir{3,200} = @(x, y, z, t, param)(0.*x.*y); 

% Pressure

% initial condition
data.u0{1}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{3}    = @(x, y, z, t, param)(0.*x.*y);

% flags
data.flag_dirichlet{1} = [5 200];
data.flag_dirichlet{2} = [5 200];
data.flag_dirichlet{3} = [5 200];

data.flag_resistance     = [2:4 6:8];
data.OutFlow_Resistance  = [1 1 1 1 1 0.5]*37500;
data.OutFlow_RefPressure = @(t) 80 * 1333;

% Model parameters
data.dynamic_viscosity   = 0.04;
data.density             = 1.06;

% Nonlinear solver
data.NonlinearSolver.tol         = 1e-6; 
data.NonlinearSolver.maxit       = 30;

% Stabilization
data.Stabilization = 'SUPG';

% Linear solver
data.LinearSolver.type           = 'MUMPS'; % MUMPS, backslash, gmres
data.LinearSolver.mumps_reordering  = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

% time 
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.005;
data.time.tf         = 1;
data.time.nonlinearity  = 'implicit';