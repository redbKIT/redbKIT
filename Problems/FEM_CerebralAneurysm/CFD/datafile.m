%% PHYSICAL DATA

% Inlet inward normal
N1  = [0 1 0];
%R1  = 1.29;%1.3558;
%XC1 = [-2.6086;   16.0572;   36.8375]; 
%Area = 5.2279;

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y);
data.force{2} = @(x, y, z, t, param)(0.*x.*y);
data.force{3} = @(x, y, z, t, param)(0.*x.*y);

% Dirichlet
%data.bcDir_t  = @(t)( 0.5*(1-cos(pi/0.5*t)).*(t<0.5) + 1 * (t>=0.5) );

data.bcDir{1,4} = @(x, y, z, t, param)( Flow_Rate( t )  * N1(1) * 1 + 0.*x.*y); 
data.bcDir{2,4} = @(x, y, z, t, param)( Flow_Rate( t )  * N1(2) * 1 + 0.*x.*y); 
data.bcDir{3,4} = @(x, y, z, t, param)( Flow_Rate( t )  * N1(3) * 1 + 0.*x.*y); 

data.bcDir{1,200} = @(x, y, z, t, param)(0.*x.*y); 
data.bcDir{2,200} = @(x, y, z, t, param)(0.*x.*y); 
data.bcDir{3,200} = @(x, y, z, t, param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y);

% initial condition
data.u0{1}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{3}    = @(x, y, z, t, param)(0.*x.*y);

% flags
data.flag_dirichlet{1} = [4 200];
data.flag_neumann{1}   = [2 3 5 6];
data.flag_pressure{1}   = [];

data.flag_dirichlet{2} = [4 200];
data.flag_neumann{2}   = [2 3 5 6];
data.flag_pressure{2}   = [];

data.flag_dirichlet{3} = [4 200];
data.flag_neumann{3}   = [2 3 5 6];
data.flag_pressure{3}   = [];

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
data.time.dt         = 0.01;
data.time.tf         = 2;
data.time.nonlinearity  = 'semi-implicit';


% % Normal Pressure
% for flag = 3 : 12
%     data.bcPrex{flag}   = @(x, y, z, t, param)( WindkesselModel( data.time.dt, t, flag, 0 ) + 0*x.*y.*z);
% end