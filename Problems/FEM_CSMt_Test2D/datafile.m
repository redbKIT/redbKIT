%DATAFILE for CSM problem

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

% Source term
data.force{1} = @(x,y,t,param)(0.*x.*y);
data.force{2} = @(x,y,t,param)(0.*x.*y);

% Dirichlet
data.bcDir{1} = @(x,y,t,param)(0.*x.*y); 
data.bcDir{2} = @(x,y,t,param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x,y,t,param)(0.*x.*y);
data.bcNeu{2} = @(x,y,t,param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x,y,t,param)(-2*(t<10) + 0.*x.*y);

% BC flag
data.flag_dirichlet{1} = [4];
data.flag_neumann{1}   = [1 2];
data.flag_pressure{1}  = [3];

data.flag_dirichlet{2} = [4];
data.flag_neumann{2}   = [1 2];
data.flag_pressure{2}  = [3];

%initial condition
data.u0{1}   = @(x,y,t,param)(0 + 0.*x.*y);
data.u0{2}   = @(x,y,t,param)(0 + 0.*x.*y);

data.du0{1}   = @(x,y,t,param)(0 + 0.*x.*y);
data.du0{2}   = @(x,y,t,param)(0 + 0.*x.*y);

% material parameters (aluminium)
data.Young   = 70000;
data.Poisson = 0.334;
data.Density = 2700;
data.model   = 'CSM';
data.Material_Model   = param{1};%'StVenantKirchhoff', 'Linear'

% Time options
data.time.t0         = 0;
data.time.dt         = 1;
data.time.tf         = param{2};
data.time.gamma      = 1/2;
data.time.beta       = 1/4;
data.time.alpha_m    = 0;
data.time.alpha_f    = 0.333;

% Linear Solver
data.LinearSolver.type              = 'backslash'; % MUMPS, backslash, gmres
data.LinearSolver.mumps_reordering  = 5;

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-6; 
data.NonLinearSolver.maxit             = 10; 

%% ROM Options
data.OfflineTraining.Solution.h5_filename       = 'Snapshots/DisplacementSnapshots.h5';
data.OfflineTraining.Solution.h5_section        = 'Displacement';
data.OfflineTraining.Solution.SamplingFrequency = 1;

data.OfflineTraining.System.h5_filename                      = 'Snapshots/SystemSnapshots.h5';
data.OfflineTraining.System.InternalForces.h5_section        = 'F_int';
data.OfflineTraining.System.InternalForces.SamplingFrequency = 1;
data.OfflineTraining.System.ExternalForces.h5_section        = 'F_ext';
data.OfflineTraining.System.ExternalForces.SamplingFrequency = 1;