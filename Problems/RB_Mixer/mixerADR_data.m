%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

a1 = 0.1;

% Source term
data.force = @(x, y, z, t, param)(0.*x.*y.*z);

% Dirichlet
data.bcDir = @(x, y, z, t, param)( 0.*x.*y.*z ...
               + param(1) *  (y>=-0.2) .*(y<1) .* (x<=1+a1/2) .*(x>=1-a1/2) .*(z>0) ... 
               + param(2) *  (y<=0.2)  .*(y>-1) .* (x<=2+a1/2) .*(x>=2-a1/2) .*(z>0)... 
               + param(3) *  (y>=-0.2) .*(y<1) .* (x<=3+a1/2) .*(x>=3-a1/2)).*(z>0);    

% Neumann
data.bcNeu = @(x, y, z, t, param)(0.*x.*y.*z);

% Robin
data.bcRob_alpha    = @(x,y,z,t,param)(0.*x);
data.bcRob_gamma    = @(x,y,z,t,param)(0.*x);
data.bcRob_fun      = @(x,y,z,t,param)(0.*x);

% BC flag
data.flag_dirichlet = [1 2];
data.flag_neumann   = [3 4];
data.flag_robin     = [];

% diffusion
data.diffusion = @(x, y, z, t, param)( 1 + 0.*x.*y);

% transport term
data.transport{1} = @(x, y, z, t, param)(  0.*x );
data.transport{2} = @(x, y, z, t, param)(  0.*x );
data.transport{3} = @(x, y, z, t, param)(  0.*x );

% reaction
data.reaction = @(x, y, z, t, param)(0 + 0.*x.*y);