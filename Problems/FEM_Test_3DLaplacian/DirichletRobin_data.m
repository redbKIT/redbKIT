%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Source term
data.force = @(x,y,t,param)(8*pi^2*(sin(2*pi*x).*sin(2*pi*y)));

% Dirichlet
data.bcDir = @(x,y,t,param)(sin(2*pi*x).*sin(2*pi*y) + 0*x.*y);  

% Neumann
data.bcNeu = @(x,y,t,param)(- 2*pi*(cos(2*pi*x).*sin(2*pi*y)) + 0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,t,param)( 0.2 + 0.*x.*y);
data.bcRob_fun      = @(x,y,t,param)( - 2*pi.*(cos(2*pi*x).*sin(2*pi*y)) + 0.2*sin(2*pi*x).*sin(2*pi*y) + 0.*x.*y);

% BC flag
data.flag_dirichlet = [1 2 3];
data.flag_neumann   = [];
data.flag_robin     = [4];

% diffusion
data.diffusion = @(x,y,t,param)(1 + 0.*x.*y);

% transport vector (first and second components)
data.trasport{1} = @(x,y,t,param)(0 + 0.*x.*y);
data.trasport{2} = @(x,y,t,param)(0 + 0.*x.*y);

% reaction
data.reaction = @(x,y,t,param)(0 + 0.*x.*y);


% exact solution
data.uexact         = @(x,y,t,param)( sin(2*pi*x).*sin(2*pi*y));
data.uxexact        = @(x,y,t,param)( 2*pi*(cos(2*pi*x).*sin(2*pi*y)));
data.uyexact        = @(x,y,t,param)( 2*pi*(sin(2*pi*x).*cos(2*pi*y)));