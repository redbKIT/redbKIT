% DATAFILE

%Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(1.*(x==-0.5)+0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1}  = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2}  = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3}  = @(x, y, z, t, param)(0.*x.*y.*z);

% bc flag
data.flag_dirichlet{1}  = [1 5 6];
data.flag_neumann{1}    = [2 3 4 7 8];

data.flag_dirichlet{2}  = [1 3 4 5 6];
data.flag_neumann{2}    = [2 7 8];

data.flag_dirichlet{3}  = [1 5 6 7 8];
data.flag_neumann{3}    = [2 3 4];

data.dynamic_viscosity = 1/150;
data.density             = 1;

% Nonlinear solver
data.NonLinearSolver.tol         = 1e-8; 
data.NonLinearSolver.maxit       = 30;

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
    
    % Preconditioner
    data.Preconditioner.type              = 'AdditiveSchwarz'; % AdditiveSchwarz, None, ILU
    data.Preconditioner.local_solver      = 'matlab_lu'; % matlab_lu, MUMPS
    data.Preconditioner.overlap_level     = 2;
    data.Preconditioner.mumps_reordering  = 7;
    data.Preconditioner.num_subdomains    = poolsize; %poolsize, number of subdomains
    
else
    % Linear Solver
    data.LinearSolver.type              = 'backslash'; % MUMPS, backslash, gmres
    data.LinearSolver.mumps_reordering  = 7;
    
    % Preconditioner
    data.Preconditioner.type              = 'None'; % AdditiveSchwarz, None, ILU
end
