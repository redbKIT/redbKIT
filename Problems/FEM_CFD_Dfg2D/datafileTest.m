% DATAFILE for DFG2
H = 0.41;
U = 1.5;
%% PHYSICAL DATA

% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.*x.*y);

% Dirichlet
data.bcDir_t  = @(t)( sin(pi*t/8) );

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

% Stabilization
data.Stabilization = param{2};

% Nonlinear solver
data.NonLinearSolver.tol         = 1e-6; 
data.NonLinearSolver.maxit       = 30;

% Linear solver
data.LinearSolver.type           = 'backslash'; % MUMPS, backslash, gmres
data.mumps_reordering            = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

% time 
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.05;
data.time.tf         = 1;
data.time.nonlinearity  = param{1};

%% Output options
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(1^2*0.1);
data.Output.DragLift.flag            = 4;
data.Output.DragLift.filename        = param{3};
