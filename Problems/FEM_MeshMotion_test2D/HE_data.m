%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Source term
data.force = @(x,y,t,param)( 0.0*x.*y);

% Dirichlet
data.bcDir = @(x,y,t,param)( rigid_motion(x,y,t,param) ) ;%deform_boundary(x,y,t,param) );  

% Neumann
data.bcNeu = @(x,y,t,param)(0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,t,param)(0.*x);
data.bcRob_gamma    = @(x,y,t,param)(0.*x);
data.bcRob_fun      = @(x,y,t,param)(0.*x);

% BC flag
data.flag_dirichlet = [1 2 3 4 5];
data.flag_neumann   = [];
data.flag_robin     = [];

% diffusion
data.diffusion = @(x,y,t,param)(1 + 0.*x.*y);

% transport vector (first and second components)
data.transport{1} = @(x,y,t,param)(0 + 0.*x.*y);
data.transport{2} = @(x,y,t,param)(0 + 0.*x.*y);

% reaction
data.reaction = @(x,y,t,param)(0 + 0.*x.*y);