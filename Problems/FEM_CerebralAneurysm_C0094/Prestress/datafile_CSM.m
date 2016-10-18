%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3} = @(x, y, z, t, param)(0 + 0.*x.*y.*z);

% Normal Pressure
data.bcPrex{210}   = @(x, y, z, t, param)(0 + 0*x.*y.*z);%dyn/cm^2 intracranial pressure produced by the cerebral spinal fluid

% BC flag
data.flag_dirichlet{1} = [1];
data.flag_neumann{1}   = [];
data.flag_pressure{1}  = [];
data.flag_FSinterface{1}  =  [200];
data.flag_robin{1}     = [210];

data.flag_dirichlet{2} = [1];
data.flag_neumann{2}   = [];
data.flag_pressure{2}  = [];
data.flag_robin{2}     = [210];
data.flag_FSinterface{2}  =  [200];

data.flag_dirichlet{3} = [1];
data.flag_neumann{3}   = [];
data.flag_pressure{3}  = [];
data.flag_robin{3}     = [210];
data.flag_FSinterface{3}  =  [200];

data.flag_dirichletNormal = [];

data.u0{1} = @(x, y, z, t, param)(0.*x.*y);
data.u0{2} = @(x, y, z, t, param)(0.*x.*y);
data.u0{3} = @(x, y, z, t, param)(0.*x.*y);

data.du0{1} = @(x, y, z, t, param)(0.*x.*y);
data.du0{2} = @(x, y, z, t, param)(0.*x.*y);
data.du0{3} = @(x, y, z, t, param)(0.*x.*y);

% material parameters
data.Material_Model = 'NeoHookean';%'StVenantKirchhoff', 'Linear','NeoHookean','RaghavanVorp'
data.Young   = 1e+7;%dyn/cm^2
data.Poisson = 0.45;
data.Density = 1.05;
data.ElasticCoefRobin = 10^3;%dyn/cm^2

data.Prestress = true;

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-6;
data.NonLinearSolver.maxit             = 10;

% Time options
data.time.t0         = 0;
data.time.dt         = 0.001; 
data.time.tf         = 1;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;
