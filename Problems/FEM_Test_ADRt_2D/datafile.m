%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Source term
data.force = @(x,y,t,param)( exp(-((x-0.5).^2 + (y-0.5).^2)./(0.07^2) ) );

% Dirichlet
data.bcDir = @(x,y,t,param)( 0.*x.*y );  

% Neumann
data.bcNeu = @(x,y,t,param)(0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,t,param)(0.*x);
data.bcRob_fun      = @(x,y,t,param)(0.*x);

% BC flag
data.flag_dirichlet = [];
data.flag_neumann   = [1 2 3 4];
data.flag_robin     = [];

% diffusion
data.diffusion = @(x,y,t,param)(0.005 + 0.*x.*y);

% transport vector (first and second components)
data.transport{1} = @(x,y,t,param)( cos(t)  + 0.*x.*y);
data.transport{2} = @(x,y,t,param)( sin(t) + 0.*x.*y);

% reaction
data.reaction = @(x,y,t,param)(0.01 + 0.*x.*y);

%data.Stabilization = 'SUPG';

% Initial condition
data.u0       = @(x,y,t,param)(0.0 + 0.*x.*y);

% time 
data.time.BDF_order = 2;
data.time.t0        = 0;
data.time.tf        = 2*pi;
data.time.dt        = 2*pi/20;