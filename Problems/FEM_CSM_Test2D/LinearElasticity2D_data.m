%DATAFILE for linear elasticity problem

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
data.bcPrex   = @(x,y,t,param)(-50 + 0.*x.*y);

% BC flag
data.flag_dirichlet{1} = [4 2];
data.flag_neumann{1}   = [1];
data.flag_pressure{1}  = [3];

data.flag_dirichlet{2} = [4 2];
data.flag_neumann{2}   = [1];
data.flag_pressure{2}  = [3];

% material parameters (aluminium)
data.Young   = 70000;
data.Poisson = 0.334;
data.model   = 'CSM';
data.Material_Model   = 'Linear';

% OutPut Options
data.Output.ComputeVonMisesStress = true;