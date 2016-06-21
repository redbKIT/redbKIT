%DATAFILE for CSM1 Test

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

% material parameters (aluminium)
data.Young   = 1.4*10^6;
data.Poisson = 0.4;
data.Density = 10^3;
data.Material_Model   = 'StVenantKirchhoff';
data.model   = 'CSM';

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-8;
data.NonLinearSolver.maxit             = 35;
