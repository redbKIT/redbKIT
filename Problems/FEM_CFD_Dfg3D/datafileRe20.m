%DATAFILE for DFG1
H = 0.41;
U = 0.45;

%% PHYSICAL DATA

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y);
data.force{2} = @(x, y, z, t, param)(0.*x.*y);
data.force{3} = @(x, y, z, t, param)(0.*x.*y);

% Dirichlet
data.bcDir_t  = @(t)( 1 );

data.bcDir{1} = @(x, y, z, t, param)( data.bcDir_t(t) * 16*U*y.*z.*(H-y).*(H-z)/H^4.*(x==0)+ 0.*x.*y); 
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y); 
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y); 


% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y);

% initial condition
data.u0{1}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{3}    = @(x, y, z, t, param)(0.*x.*y);

% flags
data.flag_dirichlet{1} = [2 10 6];
data.flag_neumann{1}   = [3];

data.flag_dirichlet{2} = [2 10 6];
data.flag_neumann{2}   = [3];

data.flag_dirichlet{3} = [2 10 6];
data.flag_neumann{3}   = [3];

% Model parameters
data.dynamic_viscosity = 1e-3;
data.density             = 1;

% Nonlinear solver
data.NonlinearSolver.tol         = 1e-6; 
data.NonlinearSolver.maxit       = 30;

% Stabilization
%data.Stabilization = 'SUPG';

% Linear solver
data.LinearSolver.type           = 'MUMPS'; % MUMPS, backslash, gmres
data.LinearSolver.mumps_reordering  = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

% time 
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.2;
data.time.tf         = 10;
data.time.nonlinearity  = 'semi-implicit';

%% Output options
% data.options.Output.computeWSS      = 0;
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(0.1*0.2^2*0.41);
data.Output.DragLift.flag            = 6;
data.Output.DragLift.filename        = 'Results/AerodynamicForces_Re20.txt';

