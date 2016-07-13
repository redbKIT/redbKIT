%DATAFILE for FSI1 - Fluid

H     = 0.41;
U_bar = 1;

% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.*x.*y);

% Dirichlet
data.bcDir_t  = @(t)( 0.5*(1-cos(pi/1*t)).*(t<1) + 1 * (t>=1));

data.bcDir{1} = @(x, y, t, param)(data.bcDir_t(t) * 1.5*U_bar*y.*(H-y)/(H/2)^2.*(x==0) + 0.*x.*y); 
data.bcDir{2} = @(x, y, t, param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

% initial condition
data.u0{1}    = @(x, y, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, t, param)(0.*x.*y);

% BC flag
data.flag_dirichlet{1}    =  [2 4 5];
data.flag_neumann{1}      =  [3];
data.flag_FSinterface{1}  =  [7 8];
data.flag_ring{1}         =  [1]; 
data.flag_ALE_fixed{1}    =  [2 4 5 3];


data.flag_dirichlet{2}    =  [2 4 5];
data.flag_neumann{2}      =  [3];
data.flag_FSinterface{2}  =  [7 8];
data.flag_ring{2}         =  [1]; 
data.flag_ALE_fixed{2}    =  [2 4 5 3];

% Model parameters
data.dynamic_viscosity  =   1;
data.density              =   10^3;
data.Stabilization        =   'SUPG';

% Nonlinear solver
data.NonLinearSolver.tol         = 1e-6; 
data.NonLinearSolver.maxit       = 30;
 
% Linear solver
data.LinearSolver.type           = 'backslash'; % MUMPS, backslash, gmres
data.mumps_reordering            = 7;


%% Time Setting
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.01; 
data.time.tf         = 10;
data.time.nonlinearity  = 'implicit';

% %% Output options
% data.options.Output.computeWSS      = 0;
% data.options.Output.computeDragLift = 1;
% data.options.Output.Drag_factor     = 2/(10^3*0.2^2*0.05);
% data.options.Output.flag_Drag       = [3];

%% Output options
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(10^3*0.2^2*0.05);
data.Output.DragLift.flag            = [5 7 8];
data.Output.DragLift.filename        = 'Results/AerodynamicForcesFSI2.txt';