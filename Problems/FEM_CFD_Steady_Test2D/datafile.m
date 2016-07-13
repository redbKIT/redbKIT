% DATAFILE

%% PHYSICAL DATA

%Source term
data.force{1} = @(x, y, t, param)(0 + 0.*x + 0.*y);
data.force{2} = @(x, y, t, param)(0 + 0.*x + 0.*y);

% Dirichlet
data.bcDir{1} = @(x, y, t, param)(4*(y-1).*(2-y).*(x==0)+0.*x.*y);    
data.bcDir{2} = @(x, y, t, param)(0 + 0.*x + 0.*y);
% Neumann
data.bcNeu{1}  = @(x, y, t, param)(0 + 0.*x + 0.*y);
data.bcNeu{2}  = @(x, y, t, param)(0 + 0.*x + 0.*y);

% bc flag
data.flag_dirichlet{1}  = [7 8];
data.flag_neumann{1}    = [9];
data.flag_dirichlet{2}  = [7 8];
data.flag_neumann{2}    = [9];

data.dynamic_viscosity   = 1/100;
data.density             = 1;

%% SOLVER OPTIONS

% Nonlinear solver
data.NonLinearSolver.first       = 'newton';
data.NonLinearSolver.switch_tol  = 1e-2; 
data.NonLinearSolver.tol         = 1e-8; 
data.NonLinearSolver.maxit       = 30;
