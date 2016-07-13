%DATAFILE for CSM3 Test

% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(-10^3*2 + 0.*x.*y);

% Dirichlet
data.bcDir{1} = @(x, y, t, param)(0.*x.*y); 
data.bcDir{2} = @(x, y, t, param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x, y, t, param)( 0.*x.*y);

% BC flag
data.flag_dirichlet{1}    =  [6];
data.flag_neumann{1}      =  [7 8];
data.flag_pressure{1}     =  [];

data.flag_dirichlet{2}    =  [6];
data.flag_neumann{2}      =  [7 8];
data.flag_pressure{2}     =  [];

data.u0{1} = @(x, y, t, param)(0.*x.*y);
data.u0{2} = @(x, y, t, param)(0.*x.*y);

data.du0{1} = @(x, y, t, param)(0.*x.*y);
data.du0{2} = @(x, y, t, param)(0.*x.*y);

% material parameters  
data.Young   = 1.4*10^6;
data.Poisson = 0.4;
data.Density = 10^3;
data.Material_Model   = 'StVenantKirchhoff';
data.model   = 'CSM';

% Linear solver (default is Matlab backslash)
data.options.LinSolver.solver            = 'backslash';%'DD_Schwarz';

% Nonlinear solver
data.options.NonlinearSolver.first       = 'newton';
data.options.NonlinearSolver.tol         = 1e-7; 
data.options.NonlinearSolver.maxIt       = 12;

% Time options
data.time.t0         = 0;
data.time.dt         = 0.005;
data.time.tf         = 1;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;
data.time.alpha_m    = 0;
data.time.alpha_f    = 0;
