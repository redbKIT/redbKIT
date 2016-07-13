%DATAFILE for DFG1
H = 0.41;
U = 0.3;
%% PHYSICAL DATA

% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.*x.*y);

% Dirichlet
data.bcDir_t  = @(t)( 1 );

data.bcDir{1} = @(x, y, t, param)( data.bcDir_t(t) * 4*U*y.*(H-y)./H^2 .* (x==0)+ 0.*x.*y); 
data.bcDir{2} = @(x, y, t, param)(0.*x.*y); 


% Neumann
data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

% initial condition
data.u0{1}    = @(x, y, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, t, param)(0.*x.*y);

% flags
data.flag_dirichlet{1} = [1 3 4];
data.flag_neumann{1}   = [2];

data.flag_dirichlet{2} = [1 3 4];
data.flag_neumann{2}   = [2];

% Model parameters
data.dynamic_viscosity = 1e-3;
data.density             = 1;

% Nonlinear solver
data.NonlinearSolver.tol         = 1e-6; 
data.NonlinearSolver.maxit       = 30;

% Stabilization
%data.Stabilization = 'SUPG';

% Linear solver
data.LinearSolver.type           = 'backslash'; % MUMPS, backslash, gmres
data.LinearSolver.mumps_reordering  = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

% time 
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.2;
data.time.tf         = 12;
data.time.nonlinearity  = 'semi-implicit';

%% Output options
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(0.2^2*0.1);
data.Output.DragLift.flag            = 4;
data.Output.DragLift.filename        = 'Results/AerodynamicForcesRe20.txt';

