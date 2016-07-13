%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Source term
data.force = @(x,y,z,t,param)( 3*pi^2* sin( pi * x ) .* sin( pi * y ) .* sin( pi * z )  );

% Dirichlet
data.bcDir = @(x,y,z,t,param)( sin( pi * x ) .* sin( pi * y ) .* sin( pi * z )); 

% Neumann
data.bcNeu = @(x,y,z,t,param)(0.*x.*y.*z);

% Robin
data.bcRob_alpha    = @(x,y,z,t,param)(0.*x);
data.bcRob_gamma    = @(x,y,z,t,param)(0.*x);
data.bcRob_fun      = @(x,y,z,t,param)(0.*x);

% BC flag
data.flag_dirichlet = [1 2 3 4 5 6];
data.flag_neumann   = [];
data.flag_robin     = [];

% diffusion
data.diffusion = @(x,y,z,t,param)(1 + 0.*x.*y);

% transport vector (first, second and third component)
data.transport{1} = @(x,y,z,t,param)(0 + 0.*x.*y);
data.transport{2} = @(x,y,z,t,param)(0 + 0.*x.*y);
data.transport{3} = @(x,y,z,t,param)(0 + 0.*x.*y);

% reaction
data.reaction = @(x,y,z,t,param)(0 + 0.*x.*y);


% exact solution
data.uexact         = @(x,y,z,t,param)( sin( pi * x ) .* sin( pi * y ) .* sin( pi * z ) );
data.uxexact        = @(x,y,z,t,param)( pi* cos( pi * x ) .* sin( pi * y ) .* sin( pi * z ) );
data.uyexact        = @(x,y,z,t,param)( pi* sin( pi * x ) .* cos( pi * y ) .* sin( pi * z ) );
data.uzexact        = @(x,y,z,t,param)( pi* sin( pi * x ) .* sin( pi * y ) .* cos( pi * z ) );

% % Linear Solver
% data.LinearSolver.type              = 'gmres'; % MUMPS, backslash, gmres
% data.LinearSolver.tol               = 1e-8; 
% data.LinearSolver.maxit             = 500; 
% data.LinearSolver.gmres_verbosity   = 5;
% data.LinearSolver.mumps_reordering  = 7;
% 
% % Preconditioner
% data.Preconditioner.type              = 'AdditiveSchwarz'; % AdditiveSchwarz, None, ILU
% data.Preconditioner.local_solver      = 'matlab_lu'; % matlab_lu, MUMPS
% data.Preconditioner.overlap_level     = 1;
% data.Preconditioner.mumps_reordering  = 7;
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     poolsize = 0;
% else
%     poolsize = poolobj.NumWorkers;
% end
% data.Preconditioner.num_subdomains    = poolsize; %poolsize, number of subdomains