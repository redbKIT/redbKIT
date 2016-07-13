%DATAFILE for FS1 - Solid

% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.0 + 0.*x.*y);

% Dirichlet
data.bcDir{1} = @(x, y, t, param)(0.*x.*y); 
data.bcDir{2} = @(x, y, t, param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x, y, t, param)(0.*x.*y);

% BC flag
data.flag_dirichlet{1}    =  [6];
data.flag_neumann{1}      =  [];
data.flag_FSinterface{1}  =  [7 8];
data.flag_pressure{1}     =  [];

data.flag_dirichlet{2}    =  [6];
data.flag_neumann{2}      =  [];
data.flag_FSinterface{2}  =  [7 8];
data.flag_pressure{2}     =  [];

data.u0{1} = @(x, y, t, param)(0.*x.*y);
data.u0{2} = @(x, y, t, param)(0.*x.*y);
data.du0{1} = @(x, y, t, param)(0.*x.*y);
data.du0{2} = @(x, y, t, param)(0.*x.*y);

% material parameters 
data.Young   = 1.4*10^6;
data.Poisson = 0.4;
data.Density = 10*10^3;
data.Material_Model   = 'StVenantKirchhoff';%'StVenantKirchhoff', Linear
data.model   = 'CSM';

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-6;
data.NonLinearSolver.maxit             = 35;

% Time options
data.time.t0         = 0;
data.time.dt         = 0.02;
data.time.tf         = 20;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;
data.time.alpha_m    = 0;
data.time.alpha_f    = 0;