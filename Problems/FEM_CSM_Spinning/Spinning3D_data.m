%DATAFILE

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

% Source term
data.force{1} = @(x,y,z,t,param)(0.*x.*y);
data.force{2} = @(x,y,z,t,param)(0.*x.*y);
data.force{3} = @(x,y,z,t,param)(-9.81 * 2700 + 0.*x.*y);

% Dirichlet
data.bcDir{1} = @(x,y,z,t,param)(0.*x.*y);
data.bcDir{2} = @(x,y,z,t,param)(0.*x.*y);
data.bcDir{3} = @(x,y,z,t,param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x,y,z,t,param)(0.*x.*y);
data.bcNeu{2} = @(x,y,z,t,param)(0.*x.*y);
data.bcNeu{3} = @(x,y,z,t,param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x,y,z,t,param)(0.*x.*y);

% BC flag
data.flag_dirichlet{1} = [];
data.flag_neumann{1}   = [1];
data.flag_pressure{1}  = [];
data.flag_clamp_points{1} = [3];

data.flag_dirichlet{2} = [];
data.flag_neumann{2}   = [1];
data.flag_pressure{2}  = [];
data.flag_clamp_points{2} = [3];

data.flag_dirichlet{3} = [];
data.flag_neumann{3}   = [1];
data.flag_pressure{3}  = [];
data.flag_clamp_points{3} = [3];

alpha = pi/3;

data.u0{1} = @(x,y,z,t,param)(0.*x.*y);
data.u0{2} = @(x,y,z,t,param)(cos(alpha)*y -sin(alpha)*z - y);
data.u0{3} = @(x,y,z,t,param)(sin(alpha)*y +cos(alpha)*z - z);

data.du0{1} = @(x,y,z,t,param)( Initial_Velocity( x, y, z, t, param, 1 ) );
data.du0{2} = @(x,y,z,t,param)( Initial_Velocity( x, y, z, t, param, 2 ) );
data.du0{3} = @(x,y,z,t,param)( Initial_Velocity( x, y, z, t, param, 3 ) );


% material parameters
data.Young   = 80*1e+9;
data.Poisson = 0.33;
data.Density = 2700;
data.Material_Model   = 'StVenantKirchhoff';
data.model   = 'CSM';

% Linear Solver
data.LinearSolver.type              = 'backslash'; % MUMPS, backslash
data.LinearSolver.mumps_reordering  = 7;

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-5; 
data.NonLinearSolver.maxit             = 15; 

% Time options
data.time.t0         = 0;
data.time.dt         = 2*1e-4;
data.time.tf         = 1;
data.time.gamma      = 0.833;
data.time.beta       = 0.444;
data.time.alpha_m    = 0;
data.time.alpha_f    = 0.333;
