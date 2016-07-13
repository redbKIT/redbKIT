%DATAFILE for FS1 - Solid

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y);
data.force{2} = @(x, y, z, t, param)(0.*x.*y);
data.force{3} = @(x, y, z, t, param)(0.*x.*y);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y);

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)(0.*x.*y);

% BC flag
data.flag_dirichlet{1}    =  [7];
data.flag_neumann{1}      =  [8];
data.flag_FSinterface{1}  =  [6];
data.flag_pressure{1}     =  [];

data.flag_dirichlet{2}    =  [7];
data.flag_neumann{2}      =  [8];
data.flag_FSinterface{2}  =  [6];
data.flag_pressure{2}     =  [];

data.flag_dirichlet{3}    =  [7 8];
data.flag_neumann{3}      =  [];
data.flag_FSinterface{3}  =  [6];
data.flag_pressure{3}     =  [];

data.u0{1} = @(x, y, z, t, param)(0.*x.*y);
data.u0{2} = @(x, y, z, t, param)(0.*x.*y);
data.u0{3} = @(x, y, z, t, param)(0.*x.*y);

data.du0{1} = @(x, y, z, t, param)(0.*x.*y);
data.du0{2} = @(x, y, z, t, param)(0.*x.*y);
data.du0{3} = @(x, y, z, t, param)(0.*x.*y);

% material parameters 
data.Young   = 2.5*10^6;
data.Poisson = 0.35;
data.Density = 0.1;
data.Material_Model   = 'NeoHookean';%'StVenantKirchhoff', Linear, NeoHookean2
data.model   = 'CSM';

mu = data.Young / (2 + 2 * data.Poisson);
lambda =  data.Young * data.Poisson /( (1+data.Poisson) * (1-2*data.Poisson) );
bulk = ( 2.0 / 3.0 ) * mu + lambda;

data.Alpha   = mu / 2;
data.Beta    = 0.0;
data.Bulk    = bulk;


% NonLinear Solver
data.NonLinearSolver.tol               = 1e-6;
data.NonLinearSolver.maxit             = 35;

% Time options
data.time.t0         = 0;
data.time.dt         = 0.005; 
data.time.tf         = 5;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;