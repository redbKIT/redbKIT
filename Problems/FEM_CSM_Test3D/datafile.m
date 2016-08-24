%DATAFILE for linear elasticity problem

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(-0.05*1 + 0.*x.*y.*z);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3} = @(x, y, z, t, param)(0 + 0.*x.*y.*z);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)(0*x.*y.*z);

% BC flag
data.flag_dirichlet{1} = [1];
data.flag_neumann{1}   = [];
data.flag_pressure{1}  = [];

data.flag_dirichlet{2} = [1];
data.flag_neumann{2}   = [];
data.flag_pressure{2}  = [];

data.flag_dirichlet{3} = [1];
data.flag_neumann{3}   = [];
data.flag_pressure{3}  = [];

% material parameters
data.Material_Model = 'StVenantKirchhoff';%'StVenantKirchhoff', 'Linear', 'NeoHookean'
data.Young   = 21.5;
data.Poisson = 0.29;
data.Density = 1;

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-6;
data.NonLinearSolver.maxit             = 25;
% data.NonLinearSolver.backtrackIter     = 3;
% data.NonLinearSolver.backtrackFactor   = 0.75;

% Solver and Preconditioner
%   If parallel pool available, use gmres with AdditiveSchwarz preconditioner
%   Otherwise use direct solver
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

if poolsize > 0
    
    % Linear Solver
    data.LinearSolver.type              = 'gmres'; % MUMPS, backslash, gmres
    data.LinearSolver.tol               = 1e-8;
    data.LinearSolver.maxit             = 500;
    data.LinearSolver.gmres_verbosity   = 5;
    data.LinearSolver.mumps_reordering  = 4;
    
    % Preconditioner
    data.Preconditioner.type         = 'AdditiveSchwarz'; % AdditiveSchwarz, None, ILU
    data.Preconditioner.local_solver = 'matlab_lu'; % matlab_lu, MUMPS
    data.Preconditioner.overlap_level     = 2;
    data.Preconditioner.mumps_reordering  = 4;
    data.Preconditioner.num_subdomains    = poolsize; %poolsize, number of subdomains
    
else
    
    % Linear Solver
    data.LinearSolver.type           = 'backslash'; % MUMPS, backslash, gmres
    
    % Preconditioner
    data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU
end

data.Output.ComputeVonMisesStress = true;