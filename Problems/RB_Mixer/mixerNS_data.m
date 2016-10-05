%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% DATAFILE for Navier-Stokes equations

%Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(z.*(0.6-z).*(y+0.2).*(0.2-y)/(0.2^2*0.3^2).*(x==-0.5)+0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1}  = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2}  = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3}  = @(x, y, z, t, param)(0.*x.*y.*z);

% bc flag
data.flag_dirichlet{1}  = [1 2];
data.flag_neumann{1}    = [3 4];

data.flag_dirichlet{2}  = [1 2];
data.flag_neumann{2}    = [3 4];

data.flag_dirichlet{3}  = [1 2 4];
data.flag_neumann{3}    = [3];

data.dynamic_viscosity   = 1/90;
data.density             = 1;

% Nonlinear solver
data.NonLinearSolver.tol     = 1e-8; 
data.NonLinearSolver.maxit   = 30;

% Linear Solver
data.LinearSolver.type    = 'backslash'; % MUMPS, backslash, gmres

% Preconditioner
data.Preconditioner.type  = 'None'; % AdditiveSchwarz, None, ILU

% Stabilization
data.Stabilization = 'SUPG';