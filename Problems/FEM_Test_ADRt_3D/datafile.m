%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Source term
data.force = @(x,y,z,t,param)( exp(-((x-0).^2 + (y-0).^2 + (z-0).^2)./(0.08^2) ) );

% Dirichlet
data.bcDir = @(x,y,z,t,param)( 0.*x.*y );  

% Neumann
data.bcNeu = @(x,y,z,t,param)(0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,z,t,param)(0.*x);
data.bcRob_gamma    = @(x,y,z,t,param)(0.*x);
data.bcRob_fun      = @(x,y,z,t,param)(0.*x);

% BC flag
data.flag_dirichlet = [];
data.flag_neumann   = [1 2 3 4 5 6];
data.flag_robin     = [];

% diffusion
data.diffusion = @(x,y,z,t,param)(0.01 + 0.*x.*y);

% transport vector (first and second components)
data.transport{1} = @(x,y,z,t,param)( cos(t)  + 0.*x.*y);
data.transport{2} = @(x,y,z,t,param)( sin(t) + 0.*x.*y);
data.transport{3} = @(x,y,z,t,param)( 0 + 0.*x.*y);

% reaction
data.reaction = @(x,y,z,t,param)(0.025 + 0.*x.*y);

%data.Stabilization = 'SUPG';

% Initial condition
data.u0       = @(x,y,z,t,param)(0.0 + 0.*x.*y);

% time 
data.time.BDF_order = 2;
data.time.t0        = 0;
data.time.tf        = 2*pi;
data.time.dt        = 2*pi/20;

% Linear Solver % MUMPS should be faster for this problem
data.LinearSolver.type              = 'backslash'; % MUMPS, backslash, gmres
data.LinearSolver.tol               = 1e-8; 
data.LinearSolver.maxit             = 500; 
data.LinearSolver.gmres_verbosity   = 1;
data.LinearSolver.mumps_reordering  = 5;

% Preconditioner
data.Preconditioner.type              = 'None'; % AdditiveSchwarz, None, ILU
data.Preconditioner.local_solver      = 'MUMPS'; % matlab_lu, MUMPS
data.Preconditioner.overlap_level     = 1;
data.Preconditioner.mumps_reordering  = 5;
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end
data.Preconditioner.num_subdomains    = 4; %poolsize, number of subdomains
data.Preconditioner.ILU_type          = 'ilutp'; 
data.Preconditioner.ILU_droptol       = 1e-3;  