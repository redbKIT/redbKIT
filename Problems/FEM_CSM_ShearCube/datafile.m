%DATAFILE for linear elasticity problem

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(0 + 0.*x.*y.*z);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3} = @(x, y, z, t, param)(1685.43 + 0.*x.*y.*z);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)(0*x.*y.*z);

% BC flag
data.flag_dirichlet{1} = [3];
data.flag_neumann{1}   = [4];
data.flag_pressure{1}  = [];

data.flag_dirichlet{2} = [3];
data.flag_neumann{2}   = [4];
data.flag_pressure{2}  = [];

data.flag_dirichlet{3} = [3];
data.flag_neumann{3}   = [4];
data.flag_pressure{3}  = [];

% material parameters 
data.Material_Model = 'StVenantKirchhoff';%'StVenantKirchhoff', 'Linear'
data.Young   = 6.6*10^4;
data.Poisson = 0.37;
data.Density = 1;

% Linear Solver
data.LinearSolver.type              = 'backslash'; % MUMPS, backslash, gmres
data.LinearSolver.mumps_reordering  = 7;

% OutPut Options
data.Output.ComputeVonMisesStress = true;