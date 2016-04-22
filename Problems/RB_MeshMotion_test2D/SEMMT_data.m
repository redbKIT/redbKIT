%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.*x.*y);

% Dirichlet
data.bcDir{1} = @(x,y,t,param)( rigid_motion(x,y,t,[1 param]) ) ;%deform_boundary(x,y,t,param) );  
data.bcDir{2} = @(x,y,t,param)( rigid_motion(x,y,t,[2 param]) ) ;%deform_boundary(x,y,t,param) );  

% Neumann
data.bcNeu{1} = @(x,y,t,param)(0.*x.*y);
data.bcNeu{2} = @(x,y,t,param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x, y, t, param)(0*x.*y);

% BC flag
data.flag_dirichlet{1} = [1 2 3 4 5];
data.flag_neumann{1}   = [];
data.flag_pressure{1}  = [];

data.flag_dirichlet{2} = [1 2 3 4 5];
data.flag_neumann{2}   = [];
data.flag_pressure{3}  = [];

data.flag_dirichlet{3} = [1 2 3 4 5];
data.flag_neumann{3}   = [];
data.flag_pressure{3}  = [];


% material parameters 
data.Material_Model = 'SEMMT';
data.Young   = 10^5;
data.Poisson = 0.4;
data.Density = 1;
data.Stiffening_power = param(3);